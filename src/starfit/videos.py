"""
Structures describing each video we are checking

Created: 1/31/25
"""
from dataclasses import dataclass, field

from starfit.good_bad_stars import goodstarsSuperTraj, goodstarsVoyagerUranus, badstarsVoyagerUranus, \
    goodstarsVoyagerNeptune, badstarsVoyagerNeptune


@dataclass(frozen=True)
class ProjectBody:
    name:str
    r:float  # Pre-encounter best estimate, used as what Blinn would use
    a:float  # Actual equatorial radius
    b:float  # Actual equatorial minor radius
    c:float  # Actual polar radius
    dt:float=0.0    # Delta-t. Adjust the post-flyby ephemeris by this many seconds.
                    # Positive moves it forward along the orbit, negative in reverse.
    parent:int=None # Orbit parent body. If set, calculate relative position between this body and
                    # parent body over a short arc of times, and treat those vectors as parent-centric vectors


solar_system_bodies:dict[int,ProjectBody]={
    10:ProjectBody(name="Sun"    ,r=0,a=0,b=0,c=0,parent=0),
     1:ProjectBody(name="Mercury",r=0,a=0,b=0,c=0,parent=0),
     2:ProjectBody(name="Venus"  ,r=0,a=0,b=0,c=0,parent=0),
     3:ProjectBody(name="Earth"  ,r=0,a=0,b=0,c=0,parent=0),
     4:ProjectBody(name="Mars"   ,r=0,a=0,b=0,c=0,parent=0),
     5:ProjectBody(name="Jupiter",r=0,a=0,b=0,c=0,parent=0),
     6:ProjectBody(name="Saturn", r=0,a=0,b=0,c=0,parent=0),
     7:ProjectBody(name="Uranus" ,r=0,a=0,b=0,c=0,parent=0),
     8:ProjectBody(name="Neptune",r=0,a=0,b=0,c=0,parent=0),
     9:ProjectBody(name="Pluto"  ,r=0,a=0,b=0,c=0,parent=0),
   }
uranus_moons:dict[int,ProjectBody]={
    799: ProjectBody(name="Uranus", r=25550, a=25559, b=25559, c=24973),
    701: ProjectBody(name="Ariel",  r=  665, a=581.1, b=577.9, c=577.7, parent=799),
    702: ProjectBody(name="Umbriel",r=  555, a=584.7, b=584.7, c=584.7, parent=799),
    703: ProjectBody(name="Titania",r=  800, a=788.9, b=788.9, c=788.9, parent=799),
    704: ProjectBody(name="Oberon", r=  815, a=761.4, b=761.4, c=761.4, parent=799),
    705: ProjectBody(name="Miranda",r=  250, a=240.4, b=234.2, c=232.9, parent=799)
}
uranus_rings=[41910,42300,42600,44800,45700,47200,47700,48300,51200] # 6,5,4,alpha,beta,eta,gamma,delta,epsilon

neptune_moons:dict[int,ProjectBody]={
        # Spice ID
        #     Name     Pre-encounter radius (km)
        #                    a        b        c (polar) from pck00010.tpc
        899: ProjectBody(name="Neptune",    r=   24781, a= 24764, b= 24764, c= 24341),
        801: ProjectBody(name="Triton",     r=  1750.0, a=1352.6, b=1352.6, c=1352.6,dt=-185.0,parent=899),
        802: ProjectBody(name="Nereid",     r=     400, a=   170, b=   170, c=   170,          parent=899),  # Radius from Neptune Travel Guide p146
        803: ProjectBody(name="Naiad N6",   r=       0, a=    29, b=    29, c=    29,          parent=899),
        804: ProjectBody(name="Thalassa N5",r=       0, a=    40, b=    40, c=    40,          parent=899),
        805: ProjectBody(name="Despina N3", r=       0, a=    74, b=    74, c=    74,          parent=899),
        806: ProjectBody(name="Galatea N4", r=       0, a=    79, b=    79, c=    79,          parent=899),
        807: ProjectBody(name="Larissa N2", r=       0, a=   104, b=   104, c=    89,          parent=899),
        808: ProjectBody(name="Proteus N1", r=       0, a=   218, b=   208, c=   201,          parent=899),
         10: ProjectBody(name="Sun"    ,    r=800000.0, a=800000, b=   208, c=   201),
        399: ProjectBody(name="Earth",      r=  6371.0, a=6371.0, b=6371.0, c=6371.0),
}
neptune_rings=[41000.0,52000.0,61000.0] # Pre-encounter estimated rings

@dataclass(frozen=True,kw_only=True)
class Project:
    framepat:str
    goodstars:set[int]=field(default_factory=set)
    badstars:set[int]=field(default_factory=set)
    bodies:dict[int,ProjectBody]=field(default_factory=dict)
    rings:list[float]=field(default_factory=list)


projects={
    "SuperTrajectory":       Project(framepat="data/frames/SuperTrajectory/frame%04d.png", goodstars=goodstarsSuperTraj,bodies=solar_system_bodies),
    "SuperTrajectoryB":      Project(framepat="data/frames/SuperTrajectoryB/frame%04d.png",goodstars=goodstarsSuperTraj,bodies=solar_system_bodies),
    "VoyagerUranusQuick":    Project(framepat='data/frames/VoyagerUranus/frame%04d.png'),
    "VoyagerUranusDetailed": Project(framepat='data/frames/VoyagerUranus/frame%04d.png',   goodstars=goodstarsVoyagerUranus,
                                     badstars=badstarsVoyagerUranus,bodies=uranus_moons,rings= uranus_rings),
    "VoyagerUranusHD":       Project(framepat='data/frames/VoyagerUranusHD/frame%04d.png', goodstars=goodstarsVoyagerUranus,
                                     badstars=badstarsVoyagerUranus, bodies=uranus_moons,rings=uranus_rings),
    "VoyagerNeptune":        Project(framepat='data/frames/VoyagerNeptune/frame%04d.png',  goodstars=goodstarsVoyagerNeptune,
                                     badstars=badstarsVoyagerNeptune,bodies=neptune_moons,rings=neptune_rings),
    "VoyagerNeptuneB":       Project(framepat='data/frames/VoyagerNeptuneB/frame%04d.png', goodstars=goodstarsVoyagerNeptune,
                                     badstars=badstarsVoyagerNeptune,bodies=neptune_moons,rings=neptune_rings)
}
