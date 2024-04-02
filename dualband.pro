PRO dualband

;; Copyright from Muhammad Aufaristama, please contact to mua2@hi.is if you one to use this script as part of your research
; Aufaristama et al 2018 (New Insights for Detecting and Deriving Thermal Properties of Lava Flow Using Infrared Satellite during 2014–2015 Effusive Eruption at Holuhraun, Iceland)

 ;;input
 OPENR, 1, 'C:\DB\btswir1'
 A1 = FLTARR(360,223) ;; This line opens the file tempertaure Band 6;; This line creates an array of size 360 (rows) by;; 223 (columns).
 OPENR, 2, 'C:\DB\bttir1'
 A2 = FLTARR(360,223) ;; This line opens the file temperature Band 10 ;; This line creates an array of size 360 (rows) by;; 223 (columns).
 OPENR, 3, 'C:\DB\TEI'
 A3 = FLTARR(360,223) ;; This line opens the file Thermal Eruption Index (TEI) ;; This line creates an array of size 360 (rows) by;; 223 (columns).
READU, 1, A1
   READU, 2, A2
      READU, 3, A3
      CLOSE, 1
   CLOSE, 2
 CLOSE, 3

  cw1 = 1.609E-6 ;; This is centralwavelength for band SWIR
  cw2 = 10.900E-6 ;; This is centralwavelength for band TIR
  BT1 = A1+273.15 ;; Conversion to Kelvin
  BT2 = A2+273.15 ;; Conversion to Kelvin
  N = A3 ;; Create TEI files into array
  Tc= (BT2-BT2)+273.15+10 ;; Cold temperature Assumption, in this case 10 degrees C and convert to Array and Kelvin
  Thstart= Tc ;; Th start same with Tc
  
  ;; Dual band method and planck inversion
  ;; Threshold from TEI
  FOR I = 0, 359 DO FOR J = 0, 222 DO BEGIN
    IF (N(i,j) GT 0.47 && N(i,j)LT 1.00) THEN Tc(i,j) = (BT2(i,j)-BT2(i,j))+273.15+85 ELSE IF (N(i,j) GT 0.2 && N(i,j)LT 0.47) THEN Tc(i,j) = (BT2(i,j)-BT2(i,j))+273.15+50 ELSE Tc(i,j) = (BT2(i,j)-BT2(i,j))+273.15+25
    Thstart(i,j) = Tc(i,j)
ENDFOR



 FOR I = 0, 359 DO FOR J = 0, 222 DO BEGIN
     IF (N(i,j) GT 0.1 && N(i,j)LT 1.00) THEN BEGIN
    IF (Thstart(i,j) LT BT1(i,j)) THEN Thstart(i,j)= BT1(i,j) ELSE Thstart(i,j) = Thstart(i,j)
    IF (Thstart(i,j) LT BT2(i,j)) THEN Thstart(i,j)= BT2(i,j) ELSE Thstart(i,j) = Thstart(i,j)

      R1 = (0.0000000000000003741*cw1^(-5.)/ (EXP(0.0143876869/(cw1*BT1))-1.));/1000000
      R2 = (0.0000000000000003741*cw2^(-5.)/ (EXP(0.0143876869/(cw2*BT2))-1.));/1000000
      Rc1 = (0.0000000000000003741*cw1^(-5.)/ (EXP(0.0143876869/(cw1*Tc))-1.));/1000000
      Rc2 = (0.0000000000000003741*cw2^(-5.)/ (EXP(0.0143876869/(cw2*Tc))-1.));/1000000
      Rhstart1 = (0.0000000000000003741*cw1^(-5.)/ (EXP(0.0143876869/(cw1*Thstart))-1.));/1000000
      Rhstart2 = (0.0000000000000003741*cw2^(-5.)/ (EXP(0.0143876869/(cw2*Thstart))-1.));/1000000
      p1 = (R1-Rc1)/(Rhstart1-Rc1)
      p2 = (R2-Rc2)/(Rhstart2-Rc2)
      Ratio = p1/p2
      
      WHILE (Ratio(i,j) NE 1.0) DO BEGIN
     IF (Ratio(i,j) GT 1.0) THEN Thstart(i,j) = Thstart(i,j) + 0.1 ELSE Thstart(i,j) = Thstart(i,j) - 0.1 
        Rh1 = (0.0000000000000003741*cw1^(-5.)/(EXP(0.0143876869/(cw1*Thstart))-1.));/1000000
        Rh2 = (0.0000000000000003741*cw2^(-5.)/(EXP(0.0143876869/(cw2*Thstart))-1.));/1000000
        p1 = (R1-Rc1)/(Rh1-Rc1)
        p2 = (R2-Rc2)/(Rh2-Rc2)
        Ratio = p1/p2
        Th1 = 0.0143876869/(cw1*alog((0.0000000000000003741* cw1^(-5.)/Rh1)+1))
        Th2 = 0.0143876869/(cw2*alog((0.0000000000000003741* cw2^(-5.)/Rh2)+1))
  

               IF (Ratio(i,j) GT 0.999) AND (Ratio(i,j) LT 1.001) THEN Ratio(i,j) = 1.0 ELSE Ratio(i,j)=Ratio(i,j)
      ENDWHILE
      Print, "Pixel Position", i,j
      Print, "Assumed Cool Component Temperature (C) = ",Tc(i,j)-273.15
      Print, "Th in Bands 1 and 2 (C) = ", Th1(i,j)-273.15,Th2(i,j)-273.15
      Print, "Solved p in Bands 1 and 2 = ", p1(i,j), p2(i,j)
      Print, "Ratio of p1/p2 = ", Ratio(i,j)


;;output
      OPENW, 4, 'C:\DB\thotnormalized3' ;; Creates a new file called thotnormalized
      WRITEU, 4, Th1-273.15 ;; Writes the data in array C into the Hot TemperatureCelcius
      CLOSE, 4 ;; This line tells IDL that we are finished reading the file.
      OPENW, 5, 'C:\DB\tcnormalized3' ;; Creates a new file called tcnormalized
      WRITEU, 5, Tc-273.15  ;; Writes the data in array C into the Cold TemperatureCelcius
      CLOSE, 5 ;; This line tells IDL that we are finished reading the file.
      OPENW, 5, 'C:\DB\fhnormalized3' ;; Creates a new file called fhnormalized
      WRITEU, 5, p1 ;; Writes the data in array C into the Pixel Fraction
      CLOSE, 5 ;; This line tells IDL that we are finished reading the file.
                
      ENDIF 
  ENDFOR

END

;; Copyright from Muhammad Aufaristama, please contact to mua2@hi.is if you one to use this script as part of your research
; Aufaristama et al 2018 (New Insights for Detecting and Deriving Thermal Properties of Lava Flow Using Infrared Satellite during 2014–2015 Effusive Eruption at Holuhraun, Iceland)
