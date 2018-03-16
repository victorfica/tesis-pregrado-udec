from features3D import features3D
from features3D.features3D_ import col_names
from prody import findPDBFiles
import os
import pandas

def make_csv(structures,outfile ,outfolder = "results_2017"):
    
    pdbs = findPDBFiles(path=structures)
    F = list(col_names().keys())
    array = []
    id_list = []
    # Si las rutas absolutas de los PDB contienen espacios, dssp no funcionara.
    for ids,path in pdbs.items():
    
        ab = features3D(pdb=path,folder_result = outfolder)
        try:
            x = ab.calc_features()
        except:
            print("Warning: {} excluded".format(ids))
            continue
        array.append(x)
        id_list.append(ids)
            
    df = pandas.DataFrame(array,id_list)      
    df.to_csv(outfile,float_format="%10.5f",header=F, index=True,sep=",")