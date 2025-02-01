#Bright star catalog, sorted by brightness.
import gzip

import numpy as np
from kwanmath.interp import linterp


def GetRA(S):
    """
    :param N: Star index in brightness order
    :return: Right Ascension in degrees
    """
    Hours=float(S[76:76+2])
    Minutes=float(S[78:78+2])
    Seconds=float(S[80:80+4])
    return (Hours+Minutes/60+Seconds/3600)*15

def GetDec(S):
    """
    :param N: Star index in brightness order
    :return: Declination in degrees
    """
    Sign   =    S[84:84+1]
    Degrees=float(S[85:85+2])
    Minutes=float(S[87:87+2])
    Seconds=float(S[89:89+2])
    return (Degrees+Minutes/60+Seconds/3600)*(-1 if Sign=="-" else 1)

def GetMag(S):
    """
    :param N: Star index in brightness order
    :return: Magnitude, lower number is brighter
    """
    return float(S[103:103+5])

def GetName(S):
    """
    :param N: Star index in brightness order
    :return: Bayer designation, Flamsteed designation, or both, followed by HR number
    """
    return S[5:5+10]+" HR"+S[1:1+4]


def GetSpectralLetter(S):
    """
    :param N: Star index in brightness order
    :return: Spectral class letter, one of OBAFGKM
    """
    return S[130:130+1]

def GetSpectralSubtype(S):
    """
    :param N: Star index in brightness order
    :return: Spectral subtype from 0 (hottest) to 9 (coolest). 10 would equivalent to 0 of the next cooler class.
    """
    return float(S[131:131+1])

def GetSpectralType(S):
    #Returns number corresponding to spectral type, O=zero, M=6
    SpectralLetter=GetSpectralLetter(S)
    if SpectralLetter=="O" or SpectralLetter=="W":
        return 0
    elif SpectralLetter=="B":
        return 1
    elif SpectralLetter=="A":
        return 2
    elif SpectralLetter=="F":
        return 3
    elif SpectralLetter=="G":
        return 4
    elif SpectralLetter=="K" or SpectralLetter=="R" or SpectralLetter=="S":
        return 5
    elif SpectralLetter=="M" or SpectralLetter=="N" or SpectralLetter=="C":
        return 6


Colors=[ #Note that these are symbolic colors, not actual colors from the blackbody curve.
    np.array((0.5 ,0.5 ,1  )), #O,W
    np.array((0.75,0.75,1  )), #B
    np.array((1   ,1   ,1  )), #A
    np.array((1   ,1   ,0.5)), #F
    np.array((1   ,1   ,0  )), #G
    np.array((1   ,0.5 ,0  )), #K,R,S
    np.array((1   ,0   ,0  )), #M,N
    np.array((0.5 ,0   ,0  ))  #(M10)
]
ColorSat=0.5
BrightMax=1
Gamma=0.4
LimitMag=6
MaxMag=-1.46

def GetColor(S):
    """
    :param N:
    :return: A 3vector for color
    """
    global LimitMag, BrightMax, Colors, BrightStarCatalog
    Mag=GetMag(S)
    #print("Mag: ",Mag)
    if(Mag>LimitMag):
        result=np.array((0,0,0)) #Star is black
    else:
        Bright=linterp(MaxMag(),BrightMax,LimitMag,0,Mag)
        Type=GetSpectralType(S)
        Subtype=GetSpectralSubtype(S)
        #print("Subtype: ",Subtype)
        #print("0 color: ",Colors[Type])
        #print("10 color: ",Colors[Type+1])
        Color=linterp(0,Colors[Type],10,Colors[Type+1],Subtype)*ColorSat
        Color=(np.array((1,1,1))*(1-ColorSat)+Color)*Bright
        result=Color**Gamma
    return result

def load_catalog():
    # Read the compressed catalog
    with gzip.open("data/starcat/catalog.gz","rt") as inf:
        lines=inf.readlines()
    # trim all the lines
    lines=[line.rstrip() for line in lines]
    # drop the lines with no data (novae etc). These ones
    # have no magnitude, so drop lines with a space in a column
    # associated with vmag. In case there are lines with
    # other valid values but not vmag, we drop them too.
    lines=[line for line in lines if line[104]!=" "]
    # Sort the catalog by brightness
    lines=list(sorted(lines,key=lambda x:float(x[103-1:107])))
    return lines

"""	
Byte-by-byte Description of file: catalog.dat
Stars are fields we care about
--------------------------------------------------------------------------------
   Bytes Format  Units     Label      Explanations
--------------------------------------------------------------------------------
   1-  4  I4     ---       HR         [1/9110]+ Harvard Revised Number
                                      = Bright Star Number
   5- 14  A10    ---       Name       Name, generally Bayer and/or Flamsteed name
  15- 25  A11    ---       DM         Durchmusterung Identification (zone in
                                      bytes 17-19)
  26- 31  I6     ---       HD         [1/225300]? Henry Draper Catalog Number
  32- 37  I6     ---       SAO        [1/258997]? SAO Catalog Number
  38- 41  I4     ---       FK5        ? FK5 star Number
      42  A1     ---       IRflag     [I] I if infrared source
      43  A1     ---       r_IRflag  *[ ':] Coded reference for infrared source
      44  A1     ---       Multiple  *[AWDIRS] Double or multiple-star code
  45- 49  A5     ---       ADS        Aitken's Double Star Catalog (ADS) designation
  50- 51  A2     ---       ADScomp    ADS number components
  52- 60  A9     ---       VarID      Variable star identification
  61- 62  I2     h         RAh1900    ?Hours RA, equinox B1900, epoch 1900.0 (1)
  63- 64  I2     min       RAm1900    ?Minutes RA, equinox B1900, epoch 1900.0 (1)
  65- 68  F4.1   s         RAs1900    ?Seconds RA, equinox B1900, epoch 1900.0 (1)
      69  A1     ---       DE-1900    ?Sign Dec, equinox B1900, epoch 1900.0 (1)
  70- 71  I2     deg       DEd1900    ?Degrees Dec, equinox B1900, epoch 1900.0 (1)
  72- 73  I2     arcmin    DEm1900    ?Minutes Dec, equinox B1900, epoch 1900.0 (1)
  74- 75  I2     arcsec    DEs1900    ?Seconds Dec, equinox B1900, epoch 1900.0 (1)
* 76- 77  I2     h         RAh        ?Hours RA, equinox J2000, epoch 2000.0 (1)
* 78- 79  I2     min       RAm        ?Minutes RA, equinox J2000, epoch 2000.0 (1)
* 80- 83  F4.1   s         RAs        ?Seconds RA, equinox J2000, epoch 2000.0 (1)
*     84  A1     ---       DE-        ?Sign Dec, equinox J2000, epoch 2000.0 (1)
* 85- 86  I2     deg       DEd        ?Degrees Dec, equinox J2000, epoch 2000.0 (1)
* 87- 88  I2     arcmin    DEm        ?Minutes Dec, equinox J2000, epoch 2000.0 (1)
* 89- 90  I2     arcsec    DEs        ?Seconds Dec, equinox J2000, epoch 2000.0 (1)
  91- 96  F6.2   deg       GLON       ?Galactic longitude (1)
  97-102  F6.2   deg       GLAT       ?Galactic latitude (1)
*103-107  F5.2   mag       Vmag       ?Visual magnitude (1)
     108  A1     ---       n_Vmag    *[ HR] Visual magnitude code
     109  A1     ---       u_Vmag     [ :?] Uncertainty flag on V
 110-114  F5.2   mag       B-V        ? B-V color in the UBV system
     115  A1     ---       u_B-V      [ :?] Uncertainty flag on B-V
 116-120  F5.2   mag       U-B        ? U-B color in the UBV system
     121  A1     ---       u_U-B      [ :?] Uncertainty flag on U-B
 122-126  F5.2   mag       R-I        ? R-I   in system specified by n_R-I
     127  A1     ---       n_R-I      [CE:?D] Code for R-I system (Cousin, Eggen)
*128-147  A20    ---       SpType     Spectral type
     148  A1     ---       n_SpType   [evt] Spectral type code
 149-154  F6.3   arcsec/yr pmRA       ?Annual proper motion in RA J2000, FK5 system
 155-160  F6.3   arcsec/yr pmDE       ?Annual proper motion in Dec J2000, FK5 system
     161  A1     ---       n_Parallax [D] D indicates a dynamical parallax,
                                      otherwise a trigonometric parallax
 162-166  F5.3   arcsec    Parallax   ? Trigonometric parallax (unless n_Parallax)
 167-170  I4     km/s      RadVel     ? Heliocentric Radial Velocity
 171-174  A4     ---       n_RadVel  *[V?SB123O ] Radial velocity comments
 175-176  A2     ---       l_RotVel   [<=> ] Rotational velocity limit characters
 177-179  I3     km/s      RotVel     ? Rotational velocity, star_vec sin i
     180  A1     ---       u_RotVel   [ :star_vec] uncertainty and variability flag on
                                      RotVel
 181-184  F4.1   mag       Dmag       ? Magnitude difference of double,
                                        or brightest multiple
 185-190  F6.1   arcsec    Sep        ? Separation of components in Dmag
                                        if occultation binary.
 191-194  A4     ---       MultID     Identifications of components in Dmag
 195-196  I2     ---       MultCnt    ? Number of components assigned to a multiple
     197  A1     ---       NoteFlag   [*] a star indicates that there is a note
                                        (file notes.dat)
--------------------------------------------------------------------------------
"""
