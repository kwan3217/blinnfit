"""
Structures describing each video we are checking

Created: 1/31/25
"""
from blinnfit.good_bad_stars import goodstarsSuperTraj, goodstarsVoyagerUranus, badstarsVoyagerUranus, \
    goodstarsVoyagerNeptune, badstarsVoyagerNeptune


projects={
    "SuperTrajectory":(675, 800,"data/frames/SuperTrajectory/frame%04d.png", goodstarsSuperTraj),
    "VoyagerUranusQuick":(201,533,'data/frames/VoyagerUranus/frame%04d.png',None,),
    "VoyagerUranusDetailed": (891, 533, 'data/frames/VoyagerUranus/frame%04d.png', goodstarsVoyagerUranus, badstarsVoyagerUranus, {
        #Spice ID
        #     Name     Pre-encounter radius (km)
        #                    a        b        c (polar) from pck00010.tpc
        799:("Uranus",25550,25559  ,25559  , 24973),
        701:("Ariel",   665,  581.1,  577.9,   577.7),
        702:("Umbriel", 555,  584.7,  584.7,   584.7),
        703:("Titania", 800,  788.9,  788.9,   788.9),
        704:("Oberon",  815,  761.4,  761.4,   761.4),
        705:("Miranda", 250,  240.4,  234.2,   232.9)
    }, None),
    "VoyagerUranusHD": (5504, 1400, 'data/frames/VoyagerUranusHD/frame%04d.png', goodstarsVoyagerUranus, badstarsVoyagerUranus, {
        #Spice ID
        #     Name     Pre-encounter radius (km)
        #                    a        b        c (polar) from pck00010.tpc
        799:("Uranus",25550,25559  ,25559  , 24973,0.0,None),
        701:("Ariel",   665,  581.1,  577.9,   577.7,0.0,799),
        702:("Umbriel", 555,  584.7,  584.7,   584.7,0.0,799),
        703:("Titania", 800,  788.9,  788.9,   788.9,0.0,799),
        704:("Oberon",  815,  761.4,  761.4,   761.4,0.0,799),
        705:("Miranda", 250,  240.4,  234.2,   232.9,0.0,799)
    }, [41910,42300,42600,44800,45700,47200,47700,48300,51200]), # 6,5,4,alpha,beta,eta,gamma,delta,epsilon
    "VoyagerNeptune": (1830, 533, 'data/frames/VoyagerNeptune/frame%04d.png', goodstarsVoyagerNeptune, badstarsVoyagerNeptune,{
        #Spice ID
        #     Name     Pre-encounter radius (km)
        #                    a        b        c (polar) from pck00010.tpc
        899:("Neptune",24781,24764  , 24764  ,24341),
        801:("Triton" , 1500, 1352.6,  1352.6,  1352.6),
        802:("Nereid" ,  400,  170  ,   170  ,   170), #Radius from Neptune Travel Guide p146
        803:("Naiad N6",   0,   29  ,    29  ,    29), #N6
        804:("Thalassa N5",0,   40  ,    40  ,    40), #N5
        805:("Despina N3" ,0,   74  ,    74  ,    74), #N3
        806:("Galatea N4" ,0,   79  ,    79  ,    79), #N4
        807:("Larissa N2" ,0,  104  ,   104  ,    89), #N2
        808:("Proteus N1" ,0,  218  ,   208  ,   201), #N1
    },[41000.0,52000.0,61000.0]
    ),
    "VoyagerNeptuneB": (
    2165, 533, 'data/frames/VoyagerNeptuneB/frame%04d.png', goodstarsVoyagerNeptune, badstarsVoyagerNeptune, {
        # Spice ID
        #     Name     Pre-encounter radius (km)
        #                    a        b        c (polar) from pck00010.tpc
        899: ("Neptune", 24781, 24764, 24764, 24341,0.0,None),
        801: ["Triton", 1750.0, 1352.6, 1352.6, 1352.6,-185.0,899],
        802: ("Nereid", 400, 170, 170, 170,0.0,899),  # Radius from Neptune Travel Guide p146
        803: ("Naiad N6", 0, 29, 29, 29,0.0,899),  # N6
        804: ("Thalassa N5", 0, 40, 40, 40,0.0,899),  # N5
        805: ("Despina N3", 0, 74, 74, 74,0.0,None),  # N3
        806: ("Galatea N4", 0, 79, 79, 79,0.0,None),  # N4
        807: ("Larissa N2", 0, 104, 104, 89,0.0,None),  # N2
        808: ("Proteus N1", 0, 218, 208, 201,0.0,None),  # N1
         10: ("Sun"    , 800000.0, 218, 208, 201,0.0,None),  # N1
        399: ("Earth",   6371.0, 218, 208, 201,0.0,None),  # N1
    }, [42800.0, 54400.0, 63400.0,68400.0]
    )
}
