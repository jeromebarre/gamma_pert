

Anthro SECTORS

agriculture
shortname: 
agr 
(1)

                agl:sector_name = "Agriculture livestock (mma)" ;
                ags:sector_name = "Agriculture soils" ;
                awb:sector_name = "Agricultural waste burning" ;

                uncertainty range D: [100%:300%]
                choose: %150
                L (correlation lenght scale) = 250km


residential
Short name:
res 
(2)
                res:sector_name = "Residential, commercial and other combustion" ;

                uncertainty range C: [50%:200%]
                choose: 75%
                L=250km

   

energy industry fugitive
Short name:
eif (3)

                ene:sector_name = "Power generation" ;
                fef:sector_name = "Fugitives" ;
                ind:sector_name = "Industrial process" ;
                slv:sector_name = "Solvents" ;

                uncertainty range B,C: [20%:60%],[50%:200%]
                choose: 50%
                L=250km 

ships
shp (4)
                shp:sector_name = "Ships" ;
        
                uncertainty range ?
                choose: 75%
                L=750km

waste
swd (5)
                swd:sector_name = "Solid waste and waste water" ;

                uncertainty range B: [20%:60%]
                choose:40%
                L=250km



transport
tra (6)
                tnr:sector_name = "Off Road transportation" ;
                tro:sector_name = "Road transportation" ;
                    
                uncertainty range C,D: [50%:200%],[100%,300%]
                choose: 100%
                L=200km


Biogenic SECTOR
bio (8)

                 megan climatology & ocean
                 uncertainty range ?
                 choose: 30%
                 L=750km

  
Fires SECTOR
fire (7)

                 GFAS
                 uncertainty range ?
                 choose: 100%
                 L=100km







