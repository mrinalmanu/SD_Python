# -*- coding: utf-8 -*-
"""
A parser for prosite data

Created on Fri Mar 15 23:00:04 2019

@author: manu
"""
from bs4 import BeautifulSoup
import urllib.request
import re

'''Url for prosite entry goes here'''
        
def search(which):
    p_id = which
    url = "https://prosite.expasy.org/{}.txt".format(p_id)
    page = urllib.request.urlopen(url)
    soup = BeautifulSoup(page, "html.parser")
        
    [s.extract() for s in soup(['style', 'script', '[document]', 'head', 'title'])]
    info = soup.getText()
            
    return(info)

    
def expand_results(info):
    psa = re.findall(r'PS.....', info)
    new_url = "https://prosite.expasy.org/cgi-bin/aligner?psa={}".format(psa[0])
    page = urllib.request.urlopen(new_url)
    soup = BeautifulSoup(page, "html.parser")
    [s.extract() for s in soup(['style', 'script', '[document]', 'head', 'title'])]
    align_info = soup.getText()
    return align_info
        
def script():
    
    print("Enter the name of protein you want to search or the prositeID")
    which = input()
    if True == bool(re.search(r"PDOC", which)):
        info = search(which)
        print(info)
    
    print("As can be see the active site of the given protein is \n \n {}".format(re.findall(r'Consensus pattern: .*\n .*', info))) 
    print("\n Do you want to see the clustal alignment as text? Yes, No")
    y = input()

    try:
        if y == 'Yes' or 'yes':
            align_info = expand_results(info)
            print(align_info)
    except:
        return "Invalid input"
        


script()







        


        


