# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 15:34:10 2022

@author: gideon.b
"""
import requests
from bs4 import BeautifulSoup
import time
import pandas
from Bio import Entrez,SeqIO
import os

Entrez.email ="gideon.bergheim@cos.uni-heidelberg.de"

#%%

def get_record_message(protein_id,out_dir):
    """
    Looks up an NCBI record and returns the information message tha tis bound 
    to the record.
    
    Parameters:
    -----------
    protein_id: string
        valid NCBI seq ID (only tested for protein Ids but others should work?)
    out_dir:    string
        folder name for output 
        
    Returns:
    --------
    message: string or Flag
        The NCBI message linke dto the id describing if the sequence is up to 
        date or has ben removed or updated.
        
    """
    if not os.path.exists("outputs/{}/{}.gb".format(out_dir,protein_id)):
        time.sleep(3)
        url="https://www.ncbi.nlm.nih.gov/protein/{}".format(protein_id)
        page = requests.get(url) 
        soup = BeautifulSoup(page.content,"html.parser")
        main = soup.find(id = "maincontent") 
        messages = main.find_all("span",class_="icon")
        if len(messages) == 0: #record unchanged
            message = "Record is unchanged."
        else:
            message = messages[0]
    else:
        message = "Record is unchanged."
    return message

def check_message(protein_id,message,out_dir,add_removed = False):
    """
    Checks the NCBI sequences message and returns and depending on its content
    returns an updated id if the record has been updated or False if the record
    has been removed or is otherwise not findable.
    returns the original id if the record is unchaged.
    
    Parameters
    ----------
    protein_id : string
        valid NCBI seq ID (only tested for protein Ids but others should work?)
    message : string or Flag
        The NCBI message linke dto the id describing if the sequence is up to 
        date or has ben removed or updated.
    out_dir : string
        folder name for output
    add_removed : bool, optional
        flag to indicate if removed records should still be listed in the 
        output excel list. The default is False.

    Returns
    -------
    out_id : string
        updated sequence id or false.

    """
    if message =="Record is unchanged.":
        out_id = protein_id
    elif  "record was removed" in message.text: 
        print("{} was removed.".format(protein_id))
        out_id = False
    elif "sequence has been updated" in message.text:
        current_id = get_current_id(message)
        print("{} was updated to: {}.".format(protein_id,current_id))
        out_id = current_id
    else:
        with open("{}/errors.log".format(out_dir),"a") as errors:
            errors.write(str(protein_id))
            errors.write("================")          
            errors.write(str(message))
            errors.write("\n\n")
        print("ERROR! {}".format(str(message)))
        out_id = False
    return out_id
    
def get_current_id(message):
    """
    If the sequence has been updated this function looks up the new sequence ID
    and returns it.
    
    Parameters
    ----------
    message : string or Flag
        The NCBI message linke dto the id describing if the sequence is up to 
        date or has ben removed or updated.
    
    Returns
    -------
    out_id : string
        updated sequence id or false.
    """
    link = message.find("a")["href"]
    current_seq_url= "https://www.ncbi.nlm.nih.gov" + link
    current_page = requests.get(current_seq_url)
    current_soup = BeautifulSoup(current_page.content,"html.parser")
    mess = current_soup.find_all("p",class_="itemid")
    for mes in mess:
        current_id = mes.text.split(": ")[1]
    return current_id
 
def retrive_sequence(protein_id,out_dir):
    """
    looks up sequence ID from NCBI and downloads it as genebank file.
    
    
    Parameters:
    -----------
    protein_id: string
        valid NCBI seq ID (only tested for protein Ids but others should work?)
    out_dir:    string
        folder name for output
    
    """
    if not os.path.exists("outputs/{}/{}.gb".format(out_dir,protein_id)):
        print("Downloading {}".format(protein_id))
        handle = Entrez.efetch(db="protein",
                                 id=protein_id,
                                 rettype="gb",
                                 retmode="text")
        with open("outputs/{}/{}.gb".format(out_dir,protein_id),"w") as f:
            f.write(handle.read())
        handle.close()
    print(protein_id,"saved as genebank.")
    
    
def check_and_download(protein_id, out_dir):
    """
    Checks a NCBI seq id to see it is up to date was removed or was updated

    Parameters:
    -----------
    protein_id: string
        valid NCBI seq ID (only tested for protein Ids but others should work?)
    out_dir:    string
        folder name for output

    Returns
    -------
    cur_id : string
        updated sequence id or false.

    """
    message = get_record_message(protein_id, out_dir)
    cur_id = check_message(protein_id, message, out_dir)
    if cur_id:
        retrive_sequence(cur_id,out_dir)
    return cur_id
        
def batch_lookup(protein_id_list, out_dir,full_report=False):
    """
    Batch look-up of NCBI sequences to check if they are up to date, were 
    removed or were updated.
    
    Parameters:
    -----------
    protein_id_list: list of strings
        list of valid NCBI ids to check
    out_dir:    string
        folder name for output    
    full_report: boolean
        flag indicating if removed records should b eincluded in the output 
        excel sheet    
    """
    
    out_dict = {}
    max_ = len(protein_id_list)
    for i,protein_id in enumerate(protein_id_list):
        print("{}/{} - Checking {}".format(i+1,max_,protein_id))
        cur_id = check_and_download(protein_id, out_dir)
        if cur_id:
            out_dict[protein_id] = cur_id
        elif full_report:
            out_dict[protein_id] = "removed"
    s = pandas.Series(out_dict)
    s.to_excel(out_dir + "\_seq_list.xlsx")
        

#%%  test case

if __name__ == "__main__":  
    test_set=["XP_001637451.2","XP_001629034.2","XP_032230999.1","WP_047149380.1"] 
    batch_lookup(test_set,"test")     


#%% checking nematostella proteins
if __name__ == "__main__": 
    ecm_list = pandas.read_csv("ecm_subset.csv",header=None)[0].tolist()
    nvec_list = [seq.id for seq in SeqIO.parse("all_sequence.fasta","fasta")]
    
    #batch_lookup(ecm_list,"ECM_subset")
    batch_lookup(nvec_list,"Nvec_proteins")
