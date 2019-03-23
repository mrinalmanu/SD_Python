# -*- coding: utf-8 -*-
"""
A parser for prosite data, also makes a fasta file to give as an input to a third party site to see conserved regions.

Created on Fri Mar 15 23:00:04 2019

@author: manu
"""
from bs4 import BeautifulSoup
import requests
import urllib.request
import re
import sys
import os


def script():
    
    print("Enter the name of protein you want to search or the prositeID e.g. PDOC00335")
    which = input()
    if True == bool(re.search(r"PDOC", which)):
        info = search_prodoc(which)
        print(info)
    elif True == bool(re.search(r"PS", which)):
        info = search_ps(which)
        print(info)
       
    print("\n Do you want to see the clustal alignment as text? Yes, No")
    y = input()

    try:
        if y == 'Yes' or 'yes':
            align_info = expand_results(info)
            print(align_info)
    except:
        return "Invalid input"
        
        '''Url for prosite entry goes here'''
        
        
def search_prodoc(which):
    p_id = which
    url = "https://prosite.expasy.org/cgi-bin/prosite/get-prodoc-entry?{}".format(p_id)
    page = urllib.request.urlopen(url)
    soup = BeautifulSoup(page, "html.parser")
        
    [s.extract() for s in soup(['style', 'script', '[document]', 'head', 'title'])]
    info = soup.getText()
            
    return(info)
    
        
def search_ps(which):
    p_id = which
    url = "https://prosite.expasy.org/cgi-bin/prosite/get-prosite-entry?{}".format(p_id)
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
    
    print("Creating fasta file")
    fasta_format = "https://prosite.expasy.org/cgi-bin/aligner?psa={}&format=FASTA".format(psa[0])
    fasta_page = urllib.request.urlopen(fasta_format)
    fasta_soup = BeautifulSoup(fasta_page, "html.parser")
    [s.extract() for s in fasta_soup(['style', 'script', '[document]', 'head', 'title'])]
    fasta_info = fasta_soup.getText()
    with open('./fasta.txt', 'a') as f1:
        f1.write(fasta_info + os.linesep)
    
    os.system('google-chrome http://www.bioinformatics.org/sms2/color_align_cons.html')
    os.system('gedit fasta.txt')
    
    return highlight(align_info)
    

def main():
    if len(sys.argv) < 1:
        print("Usage: python3 file.py")
        exit()
        
    script()
    
        
if __name__ == '__main__':
    main()
