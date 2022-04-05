#import urllib as urllib3
import re
from bs4 import *
import urllib.request
import requests
import os

def airfoiltools_database():
    aft_url = 'http://airfoiltools.com/airfoil/seligdatfile?airfoil='
    url = 'http://airfoiltools.com/search/airfoils'
    r = requests.get(url)
    temp_file = "temp_file.txt"
    with open(temp_file, 'w') as o:
        o.write(r.text)

    airfoils = []
    with open(temp_file, 'r') as inp:
        for line in inp:
            try:
                index = line.find("airfoil=")
            except:
                continue
            if index == -1:
                continue
            line2 = list(line)
            index = index + 8
            airfoil_name = ""
            while True:
                if line2[index] == "\"":
                    break
                airfoil_name += line2[index] 
                index += 1
            airfoils.append(airfoil_name)
    with open(temp_file, 'r') as inp:
        for line in inp:
            try:
                index = line.rfind("airfoil=")
            except:
                continue
            if index == -1:
                continue
            line2 = list(line)
            index = index + 8
            airfoil_name = ""
            while True:
                if line2[index] == "\"":
                    break
                airfoil_name += line2[index] 
                index += 1
            airfoils.append(airfoil_name)
    os.remove(temp_file)
    airfoils = list(dict.fromkeys(airfoils))

    i = 0
    for airfoil in airfoils:
        save_file = f"C:\\Users\\PEDRO\\Desktop\\IC DE HÉLICE\\git_control\\git_prop\\airfoils\\query_airfoiltools\\{airfoil}.dat"
        with open(save_file, 'w') as inp:
            req = requests.get(f"{aft_url}{airfoil}").text
            inp.write(req)
        i += 1
        print(f'{i} airfoil(s) saved')	

def UIUC_database():
    base_path = r"https://m-selig.ae.illinois.edu/ads/coord_seligFmt/"
    html_page = urllib.request.urlopen(base_path)																	
    soup      = BeautifulSoup(html_page,'html.parser')
    ind   = 1																																	# Initialize list of links for appending
    for link in soup.find_all('a',attrs={'href': re.compile('\.dat', re.IGNORECASE)}):				# Loop over all appropriate links on webpage
        thing = link.get('href').rsplit('/',1)[-1]
        save_path = f'airfoils\\query_UIUC\\{thing}'
        urllib.request.urlretrieve(base_path+link.get('href'), save_path)			# Get the data from the webpage, and save it to the save data file as the link name
        print("Saving file %i" %ind)																# Indicate the link that we are currently saving
        ind = ind + 1	

#UIUC_database()		
airfoiltools_database()						