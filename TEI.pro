PRO TEI

  ; Copyright from Muhammad Aufaristama, please contact to mua2@hi.is if you one to use this script as part of your research
  ; Aufaristama et al 2018 (New Insights for Detecting and Deriving Thermal Properties of Lava Flow Using Infrared Satellite during 2014–2015 Effusive Eruption at Holuhraun, Iceland)

 OPENR, 1, 'K:\Landsat 8 OLI Holuhraun\6 September 2014\IDL\toaradianceswir1'
 L1 = FLTARR(360,223)
 OPENR, 2, 'K:\Landsat 8 OLI Holuhraun\6 September 2014\IDL\toaradiancetir1'
 L2 = FLTARR(360,223)
  READU, 1, L1
  READU, 2, L2
  
A = L2*L2/((91/3)^2)
NHWLI=(L1-(L2*L2/910))/(L1+(L2*L2/(910)))*A


  OPENW, 5, 'K:\Landsat 8 OLI Holuhraun\6 September 2014\IDL\TEI'
  WRITEU, 5, TEI
  CLOSE, 5 ;; 
END