# -*- coding: utf-8 -*-
"""
Create a file containing SMILES from GZ file from PubChem
"""
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Chem.PandasTools import LoadSDF, WriteSDF, RemoveSaltsFromFrame
from contextlib import redirect_stderr
from io import StringIO
import os
from pathlib import Path
import gc
import time
#%% make Dataframe from sdf with mol object in ROMol default column, unpacking gz file first http://www.rdkit.org/docs/source/rdkit.Chem.PandasTools.html
os.chdir(r"C:\Users\Leo Liu\OneDrive - University of Toronto\Uni\CE Lab\CO2 cap")
#from PandasTools_vm import LoadSDF, WriteSDF, RemoveSaltsFromFrame #custom version of RDkit pandas tools with sanitize off
txt_path = r"D:\CE_Lab\CO2cap_Q\Filtered\sdfs\done" #input file directory
nohits_gz = []
problem_gz = []
filter1 = Chem.MolFromSmarts('[#6]1~[#6][#6](=[#8])[#6]~[#6][#6]1=[#8]') #p-benzoquinone
filter2 = Chem.MolFromSmarts('[#6]1=[#6][#6]~[#6][#6](=[#8])[#6]1=[#8]') #o-benzoquinone edge
filter3 = Chem.MolFromSmarts('c1~c~c~c~[#6](=[#8])[#6]1=[#8]') #o-benzoquinone middle
filter4 = Chem.MolFromSmarts('C1[#6](=O)~c~c~c~[#6]1(=O)') #m-quinone middle (edge meta structures are only hypothetical)
#%%
total_found = []
for files in sorted(os.listdir(txt_path)):
    if files.endswith('.sdf'):
        start_time = time.perf_counter()
        file_path = os.path.join(txt_path, files)
        print('Starting loading ', files)
        PCno = ((Path(Path(files).stem).stem))
        PubChemdf = LoadSDF(file_path, embedProps=True, 
                         removeHs=False)
        print('File loading complete ', files)
        #create folder for saving xyz files based on gz name
        folder_path = r"D:\CE_Lab\CO2cap_Q\Filtered\xyz_redone\{0}".format(PCno)
        isExist = os.path.exists(folder_path)
        if not isExist:
            os.makedirs(folder_path)
            print('folder for xyz created')
        
        # filter the molecules by substructure   

        mols_rejected1 = []
        mols_filtered1 = []

        for mol in PubChemdf['ROMol']:
            try:
                mol.UpdatePropertyCache() #use if did not sanitize
                if mol.HasSubstructMatch(filter1):
                    mols_rejected1.append(False)
                    mols_filtered1.append(mol)
                elif mol.HasSubstructMatch(filter2):
                    mols_rejected1.append(False)
                    mols_filtered1.append(mol)        
                elif mol.HasSubstructMatch(filter3):
                    mols_rejected1.append(False)
                    mols_filtered1.append(mol)
                elif mol.HasSubstructMatch(filter4):
                    mols_rejected1.append(False)
                    mols_filtered1.append(mol)         
                else:
                    mols_rejected1.append(True)
            except:
                mols_rejected1.append(True)
                continue
            
        print('Molecules passed filter:',len(mols_filtered1))
        RemoveSaltsFromFrame(PubChemdf, molCol='ROMol')

        #Remove entries containing multiple compounds
        for i, mol in enumerate(PubChemdf['ROMol']):
              frag =len(Chem.rdmolops.GetMolFrags(mol, asMols=True))
              if frag > 1:
                  mols_rejected1[i] = True
        
        PubChemdf['rejected_filter1']= mols_rejected1
        try:    
                # create 2nd dataframe with non-rejects containing original index
            PCdf_filtered_new = PubChemdf.loc[(PubChemdf['rejected_filter1'] ==False)]
            PCdf_filtered_new.reset_index(inplace=True)
            PCdf_filtered_new = PCdf_filtered_new.drop('index', axis=1)
            RemoveSaltsFromFrame(PCdf_filtered_new, molCol='ROMol')
                
                #empty out PubChemdf to free up memory
            del PubChemdf
            gc.collect()
            print('PubChemdf cleared')
    
            """ skip for problematic SDFs"""
            
            
                #GeoOpt with UFF for molecules in filtered DF, excluding any that throw an error
            print('UFF opt started')
            rdBase.LogToPythonStderr()
            mols_optimized = []
            UFFerror = []
            for mol in PCdf_filtered_new['ROMol']:
                
                try:
                    with StringIO() as buf:
                        with redirect_stderr(buf):
                            AllChem.EmbedMolecule(mol)
                            AllChem.UFFOptimizeMolecule(mol, maxIters=300)
                            res = buf.getvalue()
                            if "UFFTYPER" not in res:
                                mols_optimized.append(mol)
                            else:
                                UFFerror.append(mol)
                                mols_optimized.append(np.NaN)
                except ValueError:
                    UFFerror.append(mol)
                    mols_optimized.append(np.NaN)
                        
            PCdf_filtered_new['UFF_opt']= mols_optimized
            print('UFF errors encountered:', len(UFFerror))
                
                # final DF with UFF Value error mols removed
            PCdf_final = PCdf_filtered_new.dropna(subset =['UFF_opt'])
            PCdf_final.reset_index(inplace=True) 
            total_found.append(len(PCdf_final))
            
            #skip to here
            
            #for problematic sdfs
            PCdf_final = PCdf_filtered_new.copy()
                
                #write these mols to sdf
            WriteSDF(PCdf_final,r"D:\CE_Lab\CO2cap_Q\Filtered\sdfs\redone\{0}.sdf".format(PCno), 
                                 molColName='ROMol', properties=list(PCdf_final.columns))
            print('SDF saved')
                
                # write xyz file from the successful candidates with CID filename
            from rdkit.Chem import rdmolfiles as mf
            for i, mol in enumerate(PCdf_final['UFF_opt']):
                CID_val = PCdf_final['PUBCHEM_COMPOUND_CID'][i]
                try:
                    mf.MolToXYZFile(PCdf_final['UFF_opt'][i], r"D:\CE_Lab\CO2cap_Q\Filtered\xyz_redone\{0}\{1}.xyz".format(PCno, CID_val))
                except:
                    continue
            print('xyz saved')  
                
            del [ PCdf_filtered_new, PCdf_final ]
            gc.collect()
            print('Filtered and Final DFs cleared')
            minutes = (time.perf_counter()-start_time)/60
            print('run time in minutes: ', minutes)
        except:
            nohits_gz.append(PCno)
            pass

print('all done, molecules obtained: ', sum(total_found))
