import collections
from scipy import sparse
import numpy as np
from datetime import datetime,timedelta
import astropy.io.fits as pyfits


class CfgFile():
    def __init__(self, url):
        self.__url = url
        self.__dico = collections.OrderedDict()
        with open(self.__url, 'r') as fichier:
            content = fichier.readlines()
            fichier.close()
        content = [line.strip() for line in content]
        for line in content:
            param = line.split(" : ")
            if len(param) == 2:
                nom_param = param[0]
                value_param = param[1]
                self.__dico[nom_param] = value_param

    def __write_file(self):
        with open(self.__url, 'w') as fichier:
            for value in self.__dico.items():
                ligne = value[0] + " : " + value[1] + "\n"
                fichier.write(ligne)
            fichier.close()

    def get_param(self, nom_param):
        try:
            return self.__dico[nom_param]
        except KeyError:
            print("ERROR : Parametre " + nom_param + " non existant") 
            return ''
    def modif_param(self, nom_param, new_value):
        self.__dico[nom_param] = str(new_value)
        self.__write_file()
    
    	


def log(file,texto):
    ahora=datetime.strftime(datetime.now(),'%H:%M:%S on %d/%m/%Y')

    flog=open(file,'a')
    if ('#' in texto):
        flog.write('{0}\n'.format(texto))
    else:
        flog.write(' At {0}: {1}\n'.format(ahora,texto))
    flog.close()

def read_A_matrix(A_matrix_name):
    A_matrix = sparse.load_npz(A_matrix_name).tolil()
    [art, x_min, x_max] = A_matrix[0, -3:].toarray()[0]
    A_matrix[0, -6:] = 0

    return (A_matrix.tocsr(), art, x_min, x_max)


def header_A_matrix(A_matrix_name):
    A_matrix = sparse.load_npz(A_matrix_name).tolil()
    header = A_matrix[0, -5:-3].toarray()[0]
    return header


def date_A_matrix(A_matrix_name):
    A_matrix = sparse.load_npz(A_matrix_name).tolil()
    date_txt = str(int(A_matrix[0, -6]))
    return datetime.strptime(date_txt, "%Y%m%d%H%M")

    
def rint(nombre):
    return int(round(nombre))


def sign(x):
    return int(np.sign(x))


def checkSequence(files):
# Check if all 4 files in sequence have a True flag for data reduction
# Check if all 4 files are time ordered correctly
# Return the composite flag for DRS doing the right job

    flag=True  
    f=files[0].split('_')[1]+'_'+files[0].split('_')[2]
    t0=datetime.strptime(f,'%Y%m%d_%H%M%S')
    a=pyfits.open(files[0])
    flag=flag*a[0].header['DRS_OK']
    a.close()
    for i in range(1,4):
        f=files[i].split('_')[1]+'_'+files[i].split('_')[2]
        t=datetime.strptime(f,'%Y%m%d_%H%M%S')
        if t<t0:
            flag=False
        else:
            t0=t
        a=pyfits.open(files[i])
        flag=flag*a[0].header['DRS_OK']
        a.close()
    return  flag

def RaiseAlert(files,level,dir,texto=''):
    for f in files:
        cad=files[0].split('_')
        if len(cad)>2:
            break
    fecha0=datetime.strptime(cad[1]+'_'+cad[2],'%Y%m%d_%H%M%S')
    if  fecha0.hour<12:
        fecha0=fecha0-timedelta(1)
    fecha=datetime.strftime(fecha0,'%Y%m%d')

    alertfile=dir+'Alert_'+fecha+'.txt'
    if level==0:   #Sequence not complete
        txt='Sequence not complete: '
        for f in files:
            txt=txt+f+' '
        
    elif level==1:  #Level 3 file is not there
        txt='Level 3 file is not there: '
        txt=txt+files[0]
    elif level==2:   #No st3 at the end of  night
        txt='There  are no st3 files!!!'
    else:
        txt=texto
    txt=txt+'\n'

    if texto!='':
        txt=texto
    fout=open(alertfile,'a')
    fout.write(txt)
    fout.close()
    return 1

def cruz(img):
     
     dx=img[2050:2150,:]
     dy=img[:,2050:2150]
     img1=0.*img.copy()
     img1[0:100,:]=dx
     img1[:,0:100]=dy
     img1[100:2154,100:2148]=img[0:2054,0:2048]        #1
     img1[100:2154,2148:4196]=img[0:2054,2148:4196]    #2
     img1[2154:,100:2148]=img[2154:,0:2048]            #3
     img1[2154:,2148:]=img[2154:,2148:]                #4
     return img1
