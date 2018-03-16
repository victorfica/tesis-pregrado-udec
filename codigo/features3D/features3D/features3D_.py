# -*- coding: utf-8 -*-
from features3D import dssp_parse as ds
from collections import Counter
from glob import glob,glob1
import prody
import tempfile
import shutil
import subprocess
import os
import sys
import numpy as np
import re


AA = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
folder_result = "prediction_result"
root = os.getcwd()

def MaxASA(scale = "Tien2013"):

    """Returns the maximum accessible surface area for amino acids.
    This function returns the maximum accessible surface area (ASA)
    for the amino acids in units of square angstroms. These ASAs
    are necessary for calculating the relative solvent accessibility
    of a residue.
    There are a large variety of estimates for the exact ASA for each
    amino acid. The calling variable *scale* specifies which scale to
    use. Allowed values are the following strings:

        * 'Tien2013' : The values provided by Tien et al (Maximum
          allowed solvent accessibilities of residues in proteins),
          as defined in Table 1 in the column labeled "Theoretical"
          at http://www.plosone.org/article/info:doi/10.1371/journal.pone.0080635

    The returned variable is a dictionary *asa* keyed by each 
    upper-case one-letter amino-acid code, with the values
    being the ASA for that residue."""

    if scale == "Tien2013":

        return {"ALA":129.0,"ARG":274.0,"ASN":195.0,"ASP":193.0,"CYS":167.0,
                "GLU":223.0,"GLN":225.0,"GLY":104.0,"HIS":224.0,"ILE":197.0,
                "LEU":201.0,"LYS":236.0,"MET":224.0,"PHE":240.0,"PRO":159.0,
                "SER":155.0,"THR":172.0,"TRP":285.0,"TYR":263.0,"VAL":174.0}
    else:
        return {"ALA":106.0,"CYS":135.0,"ASP":163.0,"GLU":194.0,"PHE":197.0,
                "GLY":84.0,"HIS":184.0,"ILE":169.0,"LYS":205.0,"LEU":164.0,
                "MET":188.0,"ASN":157.0,"PRO":136.0,"GLN":198.0,"ARG":248.0,
                "SER":130.0,"THR":142.0,"VAL":142.0,"TRP":227.0,"TYR":222.0}

def make_folder(folder_result = "results"):

    if not os.path.exists(folder_result):
        os.makedirs(folder_result)
        
    return folder_result
    
def get_pdb_chain(pdbid,chain,as_object = False,folder = "results",local_path = None):

    """This function parse pdb_ids and chain name, then write the specific chain.
    Also this check for non standard aminoacid MSE and rename it as MET. To parse PDBs from a local
    folder set path in *local_path* 
    
    pdbid : a PDB identifier or a filename
    chain : The PDB chain identifier
    as_object : True for return prody object
    folder : **path** Folder to storage results
    local_path : **path** to a local Folder
    """
    
    folder_result = make_folder(folder_result = folder)
    

    if local_path is not None:
        if os.path.isdir(local_path):
            pdbs = prody.findPDBFiles(path=local_path)
            parse = prody.parsePDB(pdbs[pdbid])
        else:
            raise IOError("{0} is not a valid path".format(local_path))


    else:
        prody.pathPDBFolder(folder=folder_result)
        parse = prody.parsePDB(pdbid)
    
        
    
    protein = parse.select("protein")
    p_chain = protein.select("chain %s" %chain)
    if p_chain == None:
        
        return
        
    hv= p_chain.getHierView()
    hvc = hv[chain]
    
    
    for i,r in enumerate(hvc):
        if str(r)[:3] =="MSE":
            r.setResname("MET")
            
    if as_object:
        return hvc            
    return prody.writePDB(folder_result+"/%s_%s" %(pdbid,chain),hvc)
    
def make_dssp(pdbfile,outdir = "results"):
    """Make dssp for pdb file
    Needs modify dssp bin name 'mkdssp' in prody dssp.py"""
    
    return prody.execDSSP(pdbfile,outputdir=outdir)

def dssp_features(dsspfile):

    """ - H        Alpha helix (4-12)
        - B        Isolated beta-bridge residue
        - E        Strand
        - G        3-10 helix
        - I        pi helix
        - T        Turn
        - S        Bend
        - -        Loop """ 


    #se modifico el parser de prody
    ds_ob = ds.DSSPData()
    r = ds_ob.parseDSSP(dsspfile)

    # percentage of each ss

    ss = ds_ob.getSecStruc()

    ss_counter = Counter(ss)
    ss_total = sum(ss_counter.values())

    for s,j in ss_counter.items():
        ss_counter[s] = float("{:.5f}".format(float(j)/ss_total))
    ss_dict = dict(ss_counter)
    
    
    ###Portion of ss asa
    H_portion = 0
    B_portion = 0
    E_portion = 0
    G_portion = 0
    I_portion = 0
    T_portion = 0
    S_portion = 0
    none_portion = 0

    for x,y in zip(ss,ds_ob.acc):
        if x =='H':
            H_portion += y
        if x =='B':
            B_portion += y
        if x =='E':
            E_portion += y
        if x =='G':
            G_portion += y
        if x =='I':
            I_portion += y
        if x =='T':
            T_portion += y
        if x =='S':
            S_portion += y
        if x =='-':
            none_portion += y
    ss_portion = [H_portion,B_portion,E_portion,G_portion,I_portion,T_portion,S_portion,none_portion]
    

    features_dssp = ds_ob.getH_Pbridges()+ds_ob.getH_APbridges()+ds_ob.getH_OINJ()+ds_ob.getH_OININ1()+ds_ob.getH_OININ2()+\
                    ds_ob.getH_OININ3()+ds_ob.getH_OININ4()+ds_ob.getH_OININ5()+ds_ob.getH_OINIP0()+ds_ob.getH_OINIP1()+\
                    ds_ob.getH_OINIP2()+ds_ob.getH_OINIP3()+ds_ob.getH_OINIP4()+ds_ob.getH_OINIP5()+ds_ob.getPACC()+ss_portion
    return features_dssp

def run_pops(pdb_file,rprobe=None, coarse=True,pops='pops',residueOut=True, temp_path='/tmp/'):
    
    
    outfile = pdb_file[:-4]+".out"
    # make temp directory;
    tmp_path = tempfile.mkdtemp(dir=temp_path)

    # file name must end with '.pdb' to work with POPS
    # -> create temp file of existing pdb
    #    or write model to temp file
    handle, tmp_pdb_file = tempfile.mkstemp('.pdb', dir=tmp_path)
    os.close(handle)
    if pdb_file:
        pdb_file = os.path.abspath(pdb_file)
        shutil.copy(pdb_file, tmp_pdb_file)
    
    # chdir to temp directory, as POPS writes to current working directory
    old_dir = os.getcwd()
    try:
        os.chdir(tmp_path)

        # create the command line and run
        # catch standard out & err
        command = [pops, '--pdb',tmp_pdb_file,'--popsOut',outfile,'--noTotalOut']
        if rprobe:
            command.extend(['--rProbe', rprobe])
        if coarse:
            command.extend(['--coarse'])
        if residueOut:
            command.extend(['--residueOut'])


        p = subprocess.Popen(command, universal_newlines=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        
        try:
        # get the output, then delete the temp directory
            with open(outfile) as f:
                asa_data = f.readlines()
        except:
            #os.chdir(old_dir)
            print("Error al ejecutar POPS: Archivo {}".format(pdb_file[-10:]))
            print(err)
            return
    finally:
        os.chdir(old_dir)

    shutil.rmtree(tmp_path, ignore_errors=True)
    return asa_data


def POPS_parse(asa_data):
        
    pops_dict = {}
    start = False

    for line in asa_data:
        if (re.search("=", line) ):
            start = False

        if  (re.search('ResidNe', line) ):
            start = True
            continue

        if (start):

            res_id = int(line[17:21].strip())
            pops_dict[res_id] = {"Phob":float(line[26:32].strip()),"Phil":float(line[37:43].strip()),
                                "Total":float(line[48:54].strip()),"Resd":line[5:8].strip(),"Q(SASA)":line[59:65].strip()}
    return pops_dict

def get_rasa(pops_dict):
    # Relative accessibility (RASA)
    
    rasa = []
    MAX_ACC = MaxASA(scale ="Defecto")

    for aa in pops_dict.keys():

        resname = pops_dict[aa]["Resd"]
        asa = pops_dict[aa]["Total"]      
       
        try:
            rel_acc = asa/MAX_ACC[resname]

        except KeyError:
        # Invalid value for resname 
            rel_acc = 'NA' 
        else:

            if rel_acc > 1.0:
                rel_acc = 1.0
        rasa.append(float("{:.5f}".format(rel_acc)))
    return rasa
    
def BE_list(rasa):
    # Defining exposed and buried residues. A rasa value over 30% is a exposed residue.
    
    res_BE = []
    for rel_acc in rasa:
        
        if rel_acc*100>30.00:
            res_BE.append("E")
            #print rel_acc*100, resname,"E"
        elif rel_acc*100<30.00:
            res_BE.append("B")
            #print rel_acc*100,resname,"B"
    return res_BE
    
def Philob_rmean(pops_dict):    
    ###Phil/Phob average rasa mean
    #revisar
    MAX_ACC = MaxASA(scale ="Defecto")
    phob_rasa = []
    phil_rasa = []
    
    for aa in pops_dict.keys():
        resname = pops_dict[aa]["Resd"]
        
        if pops_dict[aa]["Phob"] > 0:

            rel_phob = pops_dict[aa]["Phob"]/MAX_ACC[resname]
            phob_rasa.append(rel_phob)
        if pops_dict[aa]["Phil"] > 0:

            rel_phil = pops_dict[aa]["Phil"]/MAX_ACC[resname]
            phil_rasa.append(rel_phil)
        
    phil_mean = sum(phil_rasa)/len(phil_rasa)
    phob_mean = sum(phob_rasa)/len(phob_rasa)
    philob_l = [float("{:.5f}".format(phil_mean)),float("{:.5f}".format(phob_mean))]

    return philob_l

def AA_rasa_mean(pops_dict,rasa):
###AA average rasa mean ###
## Revisar esta funcion, tiene errores##
    from collections import Counter

    sum_rasa = Counter()
    cntAA = Counter ()
    res_l = []
    for aa in pops_dict.keys():
        res_l.append(pops_dict[aa]["Resd"])
            
    for x,y in zip(res_l,rasa):
        sum_rasa[x] += float(y)
        cntAA[x] += 1
      
    for aa in AA:
        if aa not in list(sum_rasa):
            sum_rasa.update({aa:0})
        if aa not in list(cntAA):
            cntAA.update({aa:1}) # le puse uno para que no diera error matematico despues
            
    # Promedio rasa aa
    rasa_avg = {x: sum_rasa[x] / cntAA[x] for x in cntAA}
    
    # Entregar el mismo orden para todos los valores, los diccionarios no son ordenados. Debemos ordenar los items
    # Orden ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL
    rasa_avg = [x[1] for x in sorted(rasa_avg.items())]

    return rasa_avg

def EB_asa_(res_BE,pops_dict):
    #portion of exposed and buried asa
    
    buried_asa = 0
    exposed_asa = 0
    asa_l = []
    for aa in pops_dict.keys():
        asa_l.append(pops_dict[aa]["Total"])
            
    for x,y in zip(res_BE,asa_l):
        if x =='B':
            buried_asa += y
        elif x == 'E':
            exposed_asa += y
    EB_asa = [exposed_asa,float("{:.2f}".format(buried_asa))]
    return EB_asa

def EB_phob_(res_BE,pops_dict):
    #portion of phob Exp/Bur
    
    E_phob = 0
    B_phob = 0
    phob_l = []
    for aa in pops_dict.keys():
        phob_l.append(pops_dict[aa]["Phob"])
    for a,b in zip(res_BE,phob_l):
        if a =='E':
            E_phob += b
        elif a =='B':
            B_phob += b
    EB_phob = [float("{:.2f}".format(E_phob)),B_phob]
    return EB_phob

def EB_phil_(res_BE,pops_dict):
    #portion of phil Exp/Bur
    
    E_phil = 0
    B_phil = 0
    phil_l = []
    for aa in pops_dict.keys():
        phil_l.append(pops_dict[aa]["Phil"])
    for a,b in zip(res_BE,phil_l):
        if a =='E':
            E_phil += b
        elif a =='B':
            B_phil += b
    EB_phil = [float("{:.2f}".format(E_phil)),float("{:.2f}".format(B_phil))]
    return EB_phil

def FP_AA_EB(res_BE,pops_dict):
       ###Frequency AA in E/B
    from collections import Counter

    freqE = Counter()
    freqB = Counter()
    res_l = [pops_dict[aa]["Resd"] for aa in pops_dict.keys()]

    for a,b in zip(res_BE,res_l):
        if a =='E':
            freqE[b] += 1
        elif a =='B':
            freqB[b] += 1
    #Adding AA with 0 count 
    for aa in AA:
        if aa not in list(freqE):
            freqE.update({aa:0})
        if aa not in list(freqB):
            freqB.update({aa:0})

    freqE1 = [x[1] for x in sorted(freqE.items())]
    freqB1 = [x[1] for x in sorted(freqB.items())]

    freq_EB = freqE1 + freqB1

    ###Percentage AA in E/B


    percE = [float(x)/len(res_l) for x in freqE1]

    percB = [float(x)/len(res_l) for x in freqB1]

    perc_EB = percE + percB

    fp_EB = freq_EB + perc_EB 

    return fp_EB


def get_filepath(filetype = "pdb"):
    root = os.getcwd()
    pathresult = root+"/results/"
    #dsspfile = glob(pathresult+"*.dssp")
    filepath = glob(pathresult+"*.%s" %filetype)
    
    return filepath

def col_names():

        indx = ["n_Hpb","p_Hpb","n_Hapb","p_Hapb","nH_OINJ","pH_OINJ","nH_OIHN-1","pH_OIHN-1","nH_OIHN-2","pH_OIHN-2",
                "nH_OIHN-3","pH_OIHN-3","nH_OIHN-4","pH_OIHN-4","nH_OIHN-5","pH_OIHN-5","nH_OIHN+0","pH_OIHN+0","nH_OIHN+1",
                "pH_OIHN+1","nH_OIHN+2","pH_OIHN+2","nH_OIHN+3","pH_OIHN+3","nH_OIHN+4","pH_OIHN+4","nH_OIHN+5","pH_OIHN+5",
                "P_ASA","H_portion","B_portion","E_portion","G_portion","I_portion","T_portion","S_portion","none_portion",
                "phil_rasa","phob_rasa","A_rasa","R_rasa","N_rasa","D_rasa","C_rasa","E_rasa","Q_rasa","G_rasa","H_rasa",
                "I_rasa","L_rasa","K_rasa","M_rasa","F_rasa","P_rasa","S_rasa","T_rasa","W_rasa","Y_rasa","V_rasa","Exp_asa",
                "Bur_asa","Phob_exp","Phob_bur","Phil_exp","Phil_bur","Fq_A_E","Fq_R_E","Fq_N_E","Fq_D_E","Fq_C_E","Fq_E_E",
                "Fq_Q_E","Fq_G_E","Fq_H_E","Fq_I_E","Fq_L_E","Fq_K_E","Fq_M_E","Fq_F_E","Fq_P_E","Fq_S_E","Fq_T_E","Fq_W_E",
                "Fq_Y_E","Fq_V_E","Fq_A_B","Fq_R_B","Fq_N_B","Fq_D_B","Fq_C_B","Fq_E_B","Fq_Q_B","Fq_G_B","Fq_H_B","Fq_I_B",
                "Fq_L_B","Fq_K_B","Fq_M_B","Fq_F_B","Fq_P_B","Fq_S_B","Fq_T_B","Fq_W_B","Fq_Y_B","Fq_V_B","P_A_E","P_R_E",
                "P_N_E","P_D_E","P_C_E","P_E_E","P_Q_E","P_G_E","P_H_E","P_I_E","P_L_E","P_K_E","P_M_E","P_F_E","P_P_E",
                "P_S_E","P_T_E","P_W_E","P_Y_E","P_V_E","P_A_B","P_R_B","P_N_B","P_D_B","P_C_B","P_E_B","P_Q_B","P_G_B",
                "P_H_B","P_I_B","P_L_B","P_K_B","P_M_B","P_F_B","P_P_B","P_S_B","P_T_B","P_W_B","P_Y_B","P_V_B"]

        F = ["F%03d"% x for x in range(1,146)]
        F2indx = dict(zip(F,indx))

        return F2indx

class features3D():
    """Class for calculate structural descriptors from a PDB file. A total of 145 values are constructed
    using POPS[1] and DSSP[2] outputs.
    Method *calc_features* will Parse PDB id and Chain, write the specific chain and calculate the 
    descriptors. Also checks for non standard aminoacid MSE and renames it as MET. To parse PDBs from 
    a local folder set path in *local_path* 
    
    Parameters
    ----------

    pdbid : a PDB identifier or a filename
    chain : The PDB chain identifier
    folder_result : String. Folder name to storage results
    local_path : String. path to a local Folder
     
    """

    def __init__(self,pdb,chain = None,folder_result = "results",local_path = None):
        
        self.pdb = pdb
        self.chain = chain
        self.folder_result = folder_result
        self.local_path = local_path
        
    
    def calc_features(self):

        if self.chain is not None:
            parse_pdb = get_pdb_chain(self.pdb,self.chain,folder = self.folder_result,local_path = self.local_path)
            self.pdbfile = os.path.abspath(parse_pdb)
        else:
            folder_result = make_folder(folder_result = self.folder_result)
            self.pdbfile = os.path.abspath(self.pdb)

        dsspfile = make_dssp(self.pdbfile,outdir = self.folder_result)
        self.dsspfile = os.path.abspath(dsspfile)

        asa_data = run_pops(self.pdbfile)

        pops_dict = POPS_parse(asa_data)
        rasa = get_rasa(pops_dict)
        res_BE = BE_list(rasa)

        features = dssp_features(self.dsspfile)+Philob_rmean(pops_dict)+AA_rasa_mean(pops_dict,rasa)+EB_asa_(res_BE,pops_dict)+\
                    EB_phob_(res_BE,pops_dict)+EB_phil_(res_BE,pops_dict)+FP_AA_EB(res_BE,pops_dict)

        query_array = np.asarray(features)
        
        
        return query_array
