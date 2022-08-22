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
    if not os.path.exists("{}/{}.gb".format(out_dir,protein_id)):
        time.sleep(3)
        url="https://www.ncbi.nlm.nih.gov/protein/{}".format(protein_id)
        page = requests.get(url) 
        soup = BeautifulSoup(page.content,"html.parser")
        main = soup.find(id = "maincontent") 
        messages = main.find_all("span",class_="icon")
        if len(messages) == 0: #record unchanged
            return "Record is unchanged."
        else:
            return messages[0]
    else:
        return "Record is unchanged."

def check_message(protein_id,message,out_dir,add_removed = False):
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
     link = message.find("a")["href"]
     current_seq_url= "https://www.ncbi.nlm.nih.gov" + link
     current_page = requests.get(current_seq_url)
     current_soup = BeautifulSoup(current_page.content,"html.parser")
     mess = current_soup.find_all("p",class_="itemid")
     for mes in mess:
         current_id = mes.text.split(": ")[1]
     return current_id
 
def retrive_sequence(protein_id,out_dir):
    if not os.path.exists("{}/{}.gb".format(out_dir,protein_id)):
        print("Downloading {}".format(protein_id))
        handle = Entrez.efetch(db="protein",
                                 id=protein_id,
                                 rettype="gb",
                                 retmode="text")
        with open("{}/{}.gb".format(out_dir,protein_id),"w") as f:
            f.write(handle.read())
        handle.close()
    print(protein_id,"saved as genebank.")
    
    
def check_and_download(protein_id, out_dir):
    message = get_record_message(protein_id, out_dir)
    cur_id = check_message(protein_id, message, out_dir)
    if cur_id:
        retrive_sequence(cur_id,out_dir)
    return cur_id
        
def batch_lookup(protein_id_list, out_dir,full_report=False):
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
        

#%%    
test_set=["XP_001637451.2","XP_001629034.2","XP_032230999.1","WP_047149380.1"] 
batch_lookup(test_set,"test")     


#%%
ecm_list = pandas.read_csv("ecm_subset.csv",header=None)[0].tolist()
nvec_list = [seq.id for seq in SeqIO.parse("all_sequence.fasta","fasta")]

#batch_lookup(ecm_list,"ECM_subset")
batch_lookup(nvec_list,"Nvec_proteins")
