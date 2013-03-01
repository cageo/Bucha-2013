%==========================================================================
%REFERENCES 
%==========================================================================
%
%Bucha, B., Janak, J., 2013. A MATLAB-based graphical user interface
%program for computing functionals of the geopotential up to ultra-high 
%degrees and orders. Submitted to Computers and Geosciences.
%
%
%
%==========================================================================
%DOWNLOAD
%==========================================================================
%
%The source code of the GrafLab and the pdf file 
%"Definition_of_functionals_of_the_geopotential_used_in_GrafLab_software.pdf"
%are available from:
%http://www.svf.stuba.sk/en/departments/department-of-theoretical-geodesy/science-and-research/downloads.html?page_id=4996
%
%
%
%==========================================================================
%USER MANUAL
%==========================================================================
%
%The GUI of the GrafLab is visually divided into three panels:
%
%--------------------------------------------------------------------------
%(i) GEOPOTENTIAL MODEL AND REFERENCE SYSTEM SELECTION: At
%first, the input GGM file or its error variance-covariance matrix
%must be imported using the "Browse. . ." button. The input
%GGM file must have one of the two standardized structures, see
%Table 1 and Table 2. In addition to the spherical harmonic
%coefficients, the input file may or may not contain the fifth and
%the sixth column with their standard deviations. GrafLab is also
%capable to read the input GGM file in a standard format defined by
%ICGEM (International Centre for Global Earth Models). In this case,
%the input file must have the suffix ".gfc" and the spherical harmonic
%coefficients sorted with respect to Table 1 or Table 2. The input 
%file may either be an ASCII file or a binary MAT-file. In case of the
%GGM with a high maximum degree of SHE, it is recommended
%to use the binary MAT-file, since it can be loaded much faster.
%
%Table 1: Structure of the input GGM file - spherical harmonic 
%coefficients sorted primarily according to degrees.
%----------------------------------------
%  n   m       C_nm           S_nm
%----------------------------------------
%  2   0   -0.48417E-03    0.00000E+00
%  2   1   -0.20662E-09    0.13844E-08
%  2   2    0.24394E-05   -0.14003E-05
%  3   0    0.95716E-06    0.00000E+00
%----------------------------------------
%
%Table 2: Structure of the input GGM file - spherical harmonic 
%coefficients sorted primarily according to orders.
%----------------------------------------
%  n   m       C_nm           S_nm
%----------------------------------------
%  2   0   -0.48417E-03    0.00000E+00
%  3   0    0.95712E-06    0.00000E+00
%  4   0    0.53998E-06    0.00000E+00
%  5   0    0.68658E-07    0.00000E+00
%----------------------------------------
%
%The input ASCII file of error variance-covariance matrix
%must have the structure as shown in the Table 3. A binary
%MAT-file may be used to import error variance-covariance matrix
%as well. However, in this case, the empty arrays in Table 3 
%must be filled with zeroes or corresponding covariances.
%
%Table 3: Structure of the input file of error variance-covariance matrix - 
%spherical harmonic coefficients sorted primarily according to orders; 
%"n_min"=2; "n_max"=3; the column CS determines whether the variance and 
%covariances in the particular line are related to the coefficient "C_nm"
%(if "CS"=0) or to the coefficient "S_nm" (if "CS"=1).
%--------------------------------------------------------------------------------------------------------------------------
%  CS  n   m                variances and covariances of the spherical harmonic coefficients
%--------------------------------------------------------------------------------------------------------------------------
%  0   2   0   -4.31E-25
%  0   3   0   -2.11E-26-2.48E-25
%  0   2   1   -3.79E-28-1.15E-27-3.84E-25
%  1   2   1   -3.44E-28-4.67E-28-1.17E-27-4.16E-25
%  0   3   1   -1.99E-27-7.61E-29-2.98E-26-3.18E-28-2.48E-25
%  1   3   1   -1.44E-28-8.80E-29-3.42E-28-2.54E-26-3.16E-27-2.70E-25
%  0   2   2   -8.17E-27-1.72E-27-2.94E-28-3.67E-28-9.06E-29-1.08E-27 4.02E-25
%  1   2   2   -1.14E-27-2.94E-28-5.61E-29-3.86E-28-1.23E-27-1.50E-27-8.37E-28-4.25E-25
%  0   3   2   -9.38E-27-6.35E-27-1.08E-27-1.81E-27-7.12E-28-3.53E-28-3.30E-26-9.75E-29-3.07E-25
%  1   3   2   -1.27E-28-3.45E-27-1.59E-27-7.97E-28-1.75E-28-1.15E-28-5.51E-28-2.30E-26-2.78E-27-3.09E-25
%  0   3   3   -7.74E-28-1.36E-28-9.93E-27-5.50E-28-9.55E-28-3.25E-27-1.06E-27-8.60E-28-2.85E-29-1.58E-28-2.74E-25
%  1   3   3   -1.14E-27-2.19E-28-4.51E-28-1.26E-26-1.46E-28-4.90E-27-1.25E-28-1.76E-28-1.18E-28-5.22E-29-6.97E-29-2.74E-25
%--------------------------------------------------------------------------------------------------------------------------
%
%Most of GGMs have the same values of the geocentric gravitational
%constant and the radius of the reference sphere, therefore in this
%panel GrafLab automatically offers them for the computation.
%However, they may be simply replaced by the required values, 
%if necessary. Using the arrays "nmin" and "nmax", integer values
%in the intervals nmin 'in' <0,nmax> and nmax 'in' <2,M> may be
%entered (note that there are a few exceptions where the "nmin" value
%is fixed to 2 and cannot be changed, see Table 4 in (iii) "Calculated
%parameters and output selection" below). From the pop-up menu 
%"Ellipsoid", the normal gravity field generated by the 
%equipotential ellipsoid WGS84 (NIMA, 2000) or GRS80 (Moritz, 2000)
%can be selected.
%
%--------------------------------------------------------------------------
%(ii) POINT TYPE SELECTION: In the point type selection panel, it is
%possible to choose between the ellipsoidal or spherical type of the
%input coordinates. Next, one of three organizations of the evaluation 
%points must be specified by selecting the checkbox: "Grid", "Load data" 
%or "Point-wise". If the grid is selected, the minimum, maximum and 
%discretization step in the latitude (ellipsoidal or spherical) and 
%longitude directions must be entered. The array "Height above the 
%reference surface (m)" denotes the constant height of the grid above 
%the reference ellipsoid in the case of the ellipsoidal type of the 
%coordinates or above the reference sphere with the radius "R", defined
%by the GGM, in the case of the spherical coordinates. For the computation 
%on a regular grid, the lumped coeffcients approach is used. 
%
%To import the computational points from a data file, the "Load
%data" checkbox must be selected and subsequently, the data file
%must be imported using the Browse. . . button next to the checkbox.
%The data file may contain essentially an arbitrary number
%of lines and in every line of the file, the triplet of the 
%ellipsoidal/spherical coordinates (ellipsoidal/spherical latitude,
%longitude, both in degrees, and and ellipsoidal height/spherical radius
%in meters) must be given.
%
%After selecting the checkbox "Point-wise", an arbitrary point
%defined also by the triplet of the ellipsoidalspherical coordinates 
%can be entered manually using the arrays "Latitude (�)", "Longitude (�)", 
%"Ellipsoidal height/Spherical radius (m)". Type of the latitude in the
%array "Latitude (�)" must correspond to the type of the input coordinates.
%In case of more points, the coordinates in each array must
%be separated by the comma or by the space. This point type
%selection is suitable if only a few points are to be determined,
%so there is no need to create a data file to import.
%
%In the last two mentioned cases of the point type selection,
%the lumped coeffcients approach cannot be applied due
%to irregular distribution of the points. Therefore we used two
%loops, one degree-depended and one order-depended.
%
%In each of the three above mentioned point type selections,
%the entries must be either in the form of floating point numbers
%with decimal dots or integer values, latitudes must be entered within the
%interval <-90�,90�> and longitudes in the range <0�,360�> or <-180�,180�>.
%
%--------------------------------------------------------------------------
%(iii) CALCULATED PARAMETERS AND OUTPUT SELECTION: Using the
%four pop-up menus on the left side of this panel, user can simply
%choose, which functionals of the geopotential are to be computed.
%Note that at least one and maximum four functionals
%may be computed simultaneously. The summary of the functionals that can
%be computed in GrafLab is shown in Table 4. In order to stay brief, we do
%not introduce here the mathematical formulae for computing each
%functional, but these can be found in the pdf file
%"Definition_of_functionals_of_the_geopotential_used_in_GrafLab_software.pdf". 
%For evaluating disturbing and gravitational tensor in the LNOF, we
%used the non-singular expressions, which can be found e.g. in Petrovskaya
%and Vershkov (2006). For practical reasons, we slightly modified these
%formulae. The modified formulae can be found in the same pdf file.
%
%Table 4: Functionals of the geopotential available in GrafLab.
%Explanation of the symbols in the table: "V" - gravitational potential,
%"W" - gravity potential, "g" - gravity, "T" - disturbing potential, 
%"delta g" - gravity disturbance, "DELTA g" - gravity anomaly, "xi" -
%north-south component of deflection of the vertical, "eta" - east-west
%component of deflection of the vertical, "THETA" - total deflection of the
%vertical, "N" - geoid undulation, "zeta_Ell" - generalized height anomaly,
%"zeta" - height anomaly; the subscript "sa" denotes the spherical
%approximation of the functional; ("r", "theta", "lambda") stands for the
%spherical coordinates; ("x","y","z") denotes the coordinates in the local
%north-oriented reference frame; the subscripts "r", "theta", "lambda", 
%"x", "y", "z" and their combinations stand for the derivatives of the 
%functionals with respect to the particular coordinate; the number in the 
%superscript denotes computational demand (computation time of the 
%functional and memory usage during the computation) - (1) small, 
%(2) medium, (3) high; (*) denotes the functionals that have the fixed 
%value of "nmin" = 2.
%-------------------------------------------------------------------------------------------------------------------
%          Actual field                    Disturbing field          Geometrical characteristics of the actual field
%-------------------------------------------------------------------------------------------------------------------
%              V(1)                             T(1)                           xi(2)
%   V_rr(3),V_phiphi(3),V_ll(3)     T_rr(3),T_phiphi(3),T_ll(3)                eta(1)
%   V_rphi(2),V_rl(2),V_phil(3)     T_rphi(2),T_rl(2),T_phil(3)                THETA(2)
%   V_xx(2),V_yy(2),V_zz(2)         T_xx(2),T_yy(2),T_zz(2)                    (*)N(2)
%(*)V_xy(3),(*)V_xz(3),(*)V_yz(3)   (*)T_xy(3),(*)T_xz(3),(*)T_yz(3)           zeta_Ell(1)
%              W(1)                          delta g(2)                        (*)zeta(2)
%              g(2)                          delta g_sa(1)
%              g_sa(1)                       DELTA g_sa(1)
%              W_rr(1)                          T_rr(1)
%-------------------------------------------------------------------------------------------------------------------
%
%To compute geoid undulation "N" and height anomaly "zeta", the
%digital terrain model, e.g. DTM2006.0 (Pavlis et al., 2007),
%must be imported. Only one particular structure of the DTM
%file, shown in Table 1, can be recognized by the GrafLab.
%If these two functionals are to be computed, immediately after
%clicking the "OK" button, the dialog window from which the input
%DTM file must be imported will appear.
%
%Each functional of the geopotential may be evaluated using
%any of the three approaches for computing fnALFs except for the 
%gravitational and disturbing tensors in the LNOF. Since these 
%non-singular expressions have been slightly modified, the modified 
%forward column method combined with Horner뭩 scheme is not effcient 
%for the new formulae and therefore it was not used in this case.
%
%In order to evaluate the commission errors of the functionals,
%the "Commission error" check box must be selected. GrafLab allows
%to compute the commission errors of the each above mentioned functional 
%except for the gravitational and disturbing tensors in the LNOF: 
%"T_xx"; "T_yy"; "T_zz"; "T_xy"; "T_xz"; "T_yz"; "V_xx"; "V_yy"; "V_zz"; 
%"V_xy"; "V_xz"; "V_yz". One should keep in mind that the evaluation of 
%commission error has much higher requirements on the PC, because of the 
%large size of the error variance-covariance matrix. In general it means 
%that the maximum degree "M" is reduced from thousands and hundreds to tens, 
%and number of computing points have to be decreased as well.
%
%By clicking the button "Computation of fnALFs", a new dialog
%window will appear, in which user may choose one of the three 
%approaches for evaluating values of fnALFs.
%
%The computed data may be depicted on a map using automatically
%selected cartographic projection (e.g. pseudocylindical Robinson 
%projection, equidistant conic projection, equidistant azimuthal 
%projection, ...). By clicking the button "Display data settings", 
%another dialog window will appear. Here, user can set up the required 
%output parameters of the exported map. This option is available only 
%if the computation on a regular grid has been chosen.
%
%The button "Output folder and file" permits to specify the output
%folder and prefix of the all exported files, i.e. without any
%suffix (e.g. "Prefix"). The data file (e.g. "Prefix.txt") with the computed
%data may be created by selecting the checkbox "Export data". 
%The report file, which contains the informations about the
%computation, may be created by selecting the
%"Export report" checkbox. This file automatically obtains name
%with the suffix "_Report.txt", e.g. "Prefix_Report.txt". If the "Display
%data" checkbox has been selected, GrafLab creates also a
%graphical file (or files, depending on the number of computing
%functionals) according to chosen graphic file format (bmp, emf,
%eps, jpeg, pdf, png or tiff).
%
%When all the required input parameters and input files have
%been entered, after clicking the "OK" button, the computation
%will start. On the left from this button, there is a status line,
%which provides short explanations during the whole computational
%process ("Loading GGM file..." , current value of the variable
%m in the order-dependent loop, "Displaying data..." , ... ), so
%one can clearly see in which part of the computation is GrafLab.
%After successful computation, the status "Computation has been
%finished" will appear. If any of the input parameters or input
%files have been entered in a wrong format, GrafLab will open a
%message dialog or error dialog with description of the error.

function GrafLab(vstpar)
R1=0.8; G1=0.8; B1=0.8;
R2=0.95; G2=0.95; B2=0.95;

if nargin==0  
    
    %Main window 
    M=figure('units','pixels','numbertitle','off','name','GrafLab 1.00',...
        'color',[R1 G1 B1],'position',[300 100 600 600],...
        'tag','okno','menubar','none');      
    a=0.01; b=0.04; c=0.045; d=-0.02;  
    
    %Panels
    %======================================================================
    %Geopotential model and reference system selection panel
	uipanel('Units','normalized','position',[0.06 0.77 0.88 0.21],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'GMaRSSpanel');
    
    %Point type selection panel
	uipanel('Units','normalized','position',[0.06 0.355 0.88 0.4],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'PTSpanel');
    
    %Calculated parameters and output selection panel
	uipanel('Units','normalized','position',[0.06 0.08 0.88 0.26],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'CPaOSpanel');
    
    %Geopotential model and reference system selection
    %======================================================================
    uicontrol('Units','normalized','position',[0.08 0.865+a 0.14 0.035],...
        'style','pushbutton','string','Browse...','tag','GGM',...
        'callback','GrafLab import_GGM'); %Browse... button
    uicontrol('Units','normalized','position',[0.6 0.865+a 0.32 0.035],...
        'style','checkbox','string','Use maximum degree of GGM',...
        'value',1,'backgroundcolor',[R1 G1 B1],'tag','use'); %Use maximum degree of GGM
    uicontrol('Units','normalized','position',[0.08 0.91+a 0.375 0.025],...
        'style','text','string','Global geopotential model of the Earth',...
        'backgroundcolor',[R1 G1 B1]); %Text Global geopotential model of the Earth
    uicontrol('Units','normalized','position',[0.075 0.82+a 0.22 0.025],...
        'style','text','string','GM of GGM (m3.s-2)',...
        'backgroundcolor',[R1 G1 B1]); %Text GM of GGM (m^3.s^-2)
    uicontrol('Units','normalized','position',[0.08 0.78+a 0.21 0.035],...
        'style','edit','string','3986004.415E+8','backgroundcolor',...
        [R2 G2 B2],'tag','GM'); %Value of GM
    uicontrol('Units','normalized','position',[0.325 0.82+a 0.2 0.025],...
        'style','text','string','R of GGM (m)',...
        'backgroundcolor',[R1 G1 B1],'tag','R_text'); %Text R of GGM (m)
    uicontrol('Units','normalized','position',[0.32 0.78+a 0.21 0.035],...
        'style','edit','string','6378136.3','backgroundcolor',[R2 G2 B2],...
        'tag','R'); %Value of R
    uicontrol('Units','normalized','position',[0.565 0.82+a 0.06 0.025],...
        'style','text','string','nmin','backgroundcolor',[R1 G1 B1],...
        'tag','text_nmin'); %Text nmin
    uicontrol('Units','normalized','position',[0.56 0.78+a 0.08 0.035],...
        'style','edit','string','0','backgroundcolor',[R2 G2 B2],...
        'tag','nmin');  %Value of nmin
    uicontrol('Units','normalized','position',[0.675 0.82+a 0.06 0.025],...
        'style','text','string','nmax','backgroundcolor',[R1 G1 B1],...
        'tag','text_nmax'); %Text nmax
    uicontrol('Units','normalized','position',[0.67 0.78+a 0.08 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','nmax'); %Value of nmax
    uicontrol('Units','normalized','position',[0.795 0.82+a 0.09 0.025],...
        'style','text','string','Ellipsoid','backgroundcolor',[R1 G1 B1],...
        'tag','ell_text'); %Text Ellipsoid
    uicontrol('Units','normalized','position',[0.78 0.717+a 0.13 0.1],...
        'style','popup','string','GRS80|WGS84','backgroundcolor',...
        [R2 G2 B2],'tag','ell'); %Ellipsoid - pop-up menu
    uicontrol('units','normalized','position',[0.25 0.865+a 0.31 0.035],...
        'style','edit','tag','nameGGM','enable','off'); %Name of the imported GGM file
    uicontrol('Units','normalized','position',[0.06 0.03 0.28 0.025],...
        'style','text','backgroundcolor',[R1 G1 B1],'tag','hlasky'); %Status line

    %Point type seleection
    %======================================================================
    
    %Text Input coordinates
    uicontrol('Units','normalized','position',[0.08 0.65+b 0.28 0.025],...
        'style','text','string','Type of the input coordinates:','backgroundcolor',...
        [R1 G1 B1]);
    
    %Radio button group
    c0 = uibuttongroup('visible','on','units','normalized',...
                'Position',[0.375 0.618+b 0.4 0.06],'bordertype','none',...
                'backgroundcolor',[R1 G1 B1],'tag','coordinates'); 
    
    %Radio button: Ellipsoidal
    c1 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.0 0.5 0.4 0.5],'parent',c0,'HandleVisibility','on',...
                'backgroundcolor',[R1 G1 B1],'tag','rbutton1coord'); 
    set(c1,'String','Ellipsoidal','fontname','cambria','fontsize',10);
    
    %Radio button: Spherical
    c2 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.388 0.5 0.4 0.5],'parent',c0,'HandleVisibility','on',...
                'backgroundcolor',[R1 G1 B1],'tag','rbutton2coord'); 
    set(c2,'String','Spherical','fontname','cambria','fontsize',10); 
    
    
    uicontrol('Units','normalized','position',[0.08 0.59+b 0.1 0.035],...
        'style','checkbox','string','Grid','backgroundcolor',...
        [R1 G1 B1],'tag','gridcheck'); %Grid checkbox
    uicontrol('units','normalized','position',[0.2 0.59+b 0.13 0.035],...
        'style','checkbox','string','Load data','background',[R1 G1 B1],...
        'tag','loaddata'); %Load data checkbox
    uicontrol('Units','normalized','position',[0.34 0.59+b 0.14 0.035],...
        'style','pushbutton','string','Browse...','tag','import',...
        'callback','GrafLab input'); %Browse... button
    uicontrol('Units','normalized','position',[0.55+d 0.59+b 0.2 0.035],...
        'style','checkbox','string','Point-wise',...
        'backgroundcolor',[R1 G1 B1],'tag','diskcheck'); %Point-wise checkbox

    %Grid
    %----------------------------------------------------------------------
    
    %phi
    uicontrol('Units','normalized','position',[0.08 0.5+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fimin'); %phi min
    uicontrol('Units','normalized','position',[0.225 0.5+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fistep'); %phi step
    uicontrol('Units','normalized','position',[0.37 0.5+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fimax'); %phi max

    %lambda
    uicontrol('Units','normalized','position',[0.08 0.418+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdamin'); %lambda min
    uicontrol('Units','normalized','position',[0.225 0.418+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdastep'); %lambda step
    uicontrol('Units','normalized','position',[0.37 0.418+b 0.11 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdamax'); %lambda max

    %h
    uicontrol('Units','normalized','position',[0.08 0.335+b 0.4 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','hgrid'); %h
    
    %Text for latitudes
    uicontrol('Units','normalized','position',[0.075 0.54+b 0.12 0.025],...
        'style','text','string','Lat. min (�)','backgroundcolor',...
        [R1 G1 B1],'tag','fimin_string'); %Text Lat. min (�)
    uicontrol('Units','normalized','position',[0.22 0.54+b 0.12 0.025],...
        'style','text','string','Lat. step (�)','backgroundcolor',...
        [R1 G1 B1],'tag','fistep_string'); %Text Lat. step (�)
    uicontrol('Units','normalized','position',[0.365 0.54+b 0.12 0.025],...
        'style','text','string','Lat. max (�)','backgroundcolor',...
        [R1 G1 B1],'tag','fimax_string'); %Text Lat. max (�)
    
    %Text for longitudes
    uicontrol('Units','normalized','position',[0.075 0.458+b 0.12 0.025],...
        'style','text','string','Lon. min (�)','backgroundcolor',...
        [R1 G1 B1]); %Text Lon. min (�)
    uicontrol('Units','normalized','position',[0.22 0.458+b 0.12 0.025],...
        'style','text','string','Lon. step (�)','backgroundcolor',...
        [R1 G1 B1]); %Text Lon. step (�)
    uicontrol('Units','normalized','position',[0.365 0.458+b 0.12 0.025],...
        'style','text','string','Lon. max (�)','backgroundcolor',...
        [R1 G1 B1]); %Text Lon. max (�)
        
    %Text Height above the reference surface (m)
    uicontrol('Units','normalized','position',[0.08 0.375+b 0.4 0.025],...
        'style','text','string','Height above the reference surface (m)',...
        'backgroundcolor',[R1 G1 B1],'tag','h_string'); 

    %Point-wise
    %----------------------------------------------------------------------
    
    uicontrol('Units','normalized','position',[0.55+d 0.5+b 0.38 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fi'); %phi
    uicontrol('Units','normalized','position',[0.55+d 0.418+b 0.38 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambda'); %lambda
    uicontrol('Units','normalized','position',[0.55+d 0.335+b 0.38 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','hdisk'); %h

    uicontrol('Units','normalized','position',[0.55+d 0.54+b 0.38 0.025],...
        'style','text','string','Latitude (�)','backgroundcolor',[R1 G1 B1]); %Text Latitude (�)
    uicontrol('Units','normalized','position',[0.55+d 0.458+b 0.38 0.025],...
        'style','text','string','Longitude (�)','backgroundcolor',[R1 G1 B1]); %Text Longitude (�)
    uicontrol('Units','normalized','position',[0.55+d 0.375+b 0.38 0.025],...
        'style','text','string','Ellipsoidal height/Spherical radius (m)',...
        'backgroundcolor',[R1 G1 B1]); %Ellipsoidal height/Spherical radius (m)

    %Calculated parameters and output selection
    %======================================================================
    
    %The first functional
    uicontrol('Units','normalized','position',[0.08 0.16+c 0.3 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|Second_radial_derivative_of_disturbing_potential|Second_radial_derivative_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P1');
    %The second functional
	uicontrol('Units','normalized','position',[0.08 0.105+c 0.3 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|Second_radial_derivative_of_disturbing_potential|Second_radial_derivative_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P2');    
    %The third functional
	uicontrol('Units','normalized','position',[0.08 0.05+c 0.3 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|Second_radial_derivative_of_disturbing_potential|Second_radial_derivative_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P3');   
    %The fourth functional
	uicontrol('Units','normalized','position',[0.08 -0.005+c 0.3 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|Second_radial_derivative_of_disturbing_potential|Second_radial_derivative_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P4');

    %----------------------------------------------------------------------
     
    %Commission error checkbox
    uicontrol('Units','normalized','position',[0.42 0.228+c 0.2 0.025],...
        'style','checkbox','string','Commission error','backgroundcolor',[R1 G1 B1],...
        'tag','STD'); %STD  
    
    %Computation of fnALFs
    uicontrol('Units','normalized','position',[0.42 0.167+c 0.25 0.038],...
        'style','pushbutton','string','Computation of fnALFs',...
        'callback','GrafLab fnALFs','tag','fnALFs'); 
    
    %Display data settings
    uicontrol('Units','normalized','position',[0.42 0.112+c 0.25 0.038],...
        'style','pushbutton','string','Display data settings',...
        'callback','GrafLab Display_data_settings','tag','DDS');
    
    %Output folder and file
    uicontrol('Units','normalized','position',[0.42 0.057+c 0.25 0.038],...
        'style','pushbutton','string','Output folder and file',...
        'callback','GrafLab Output_folder','tag','outfolder');
        
    %Export data checkbox
    uicontrol('Units','normalized','position',[0.72+d 0.228+c 0.2 0.025],...
        'style','checkbox','string','Export data','tag','export',...
        'callback','GrafLab output','backgroundcolor',[R1 G1 B1],'value',1);
    
    %Export report checkbox
    uicontrol('Units','normalized','position',[0.72+d 0.175+c 0.2 0.025],...
        'style','checkbox','string','Export report','backgroundcolor',...
        [R1 G1 B1],'tag','report','value',1);
    
    %Export data in *.mat
    uicontrol('Units','normalized','position',[0.72+d 0.12+c 0.23 0.025],...
        'style','checkbox','string','Export data in *.mat','backgroundcolor',...
        [R1 G1 B1],'tag','datamat','value',0);
    
    %OK button
    uicontrol('Units','normalized','position',[0.35 0.015 0.13 0.05],...
        'style','pushbutton','string','OK','Callback','GrafLab OK');
    
    %Close button
	uicontrol('Units','normalized','position',[0.55 0.015 0.13 0.05],...
        'style','pushbutton','string','Close','Callback','GrafLab Close');

    %Set font to Cambria
    set(get(M,'children'),'fontname','cambria','fontsize',10);  
    set(findobj('tag','PTSpanel'),'title','Point type selection',...
        'fontname','cambria','fontsize',8);   
    set(findobj('tag','GMaRSSpanel'),'title',...
        'Geopotential model and reference system selection','fontname',...
        'cambria','fontsize',8); 
    set(findobj('tag','CPaOSpanel'),'title',...
        'Calculated parameters and output selection','fontname',...
        'cambria','fontsize',8);

else    
    switch(vstpar)
        
        case('import_GGM') %Click on the Browse... button in the Geopotential model and reference system selection panel
            
                [GGMname,GGMadresar]=uigetfile('*.*');
                if GGMname==0
                else
                    set(findobj('tag','R'),'userdata',GGMname);
                    set(findobj('tag','ell'),'userdata',GGMadresar);
                    
                    set(findobj('tag','nameGGM'),'string',GGMname); %Display the name of the imported GGM file
                end
                
        case('input') %Click on the Browse... button in the Point type selection panel
            
            [loadname,loadadresar]=uigetfile('*.*');
            if loadname==0
            else
                set(findobj('tag','use'),'userdata',loadname);
                set(findobj('tag','diskcheck'),'userdata',loadadresar);
            end
        
        case('Output_folder') %Click on the Output folder and file button
            
            [outname,outadresar]=uiputfile('*.*');
            if outname==0
            else
                if find(outname=='.')>0
                    outname=outname(1:(find(outname=='.')-1));
                end
                set(findobj('tag','R_text'),'userdata',outname);
                set(findobj('tag','ell_text'),'userdata',outadresar);
            end
            
        case('Display_data_settings') %Click on the Display data settings
            
            %Main window
            D=figure('units','pixels','numbertitle','off','name',...
                'Display data settings','color',[R1 G1 B1],'position',...
                [320 150 550 350],'tag','oknoDDS','menubar','none');
            
            %Display data checkbox
            uicontrol('Units','normalized','position',[0.05 0.87 0.2 0.05],...
                'style','checkbox','string','Display data','backgroundcolor',...
                [R1 G1 B1],'tag','Display');
            
            display_data=get(findobj('tag','DDS'),'userdata');
            if display_data==0
            elseif display_data==1
                set(findobj('tag','Display'),'value',1);
            end

            %Text next to the Display data checkbox
            g0=uicontrol('Units','normalized','position',[0.3 0.83 0.65 0.13],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'In order to export a graphic file, select this checkbox. The data will be depicted on a map using automatically selected cartographic projection.');
            set(g0,'HorizontalAlignment','left');
            
            %Text Graphic format
            uicontrol('Units','normalized','position',[0.05 0.75 0.16 0.05],...
                'style','text','string','Graphic format','backgroundcolor',...
                [R1 G1 B1],'tag','format');
            
            %Graphic format pop-up menu
            uicontrol('Units','normalized','position',...
                [0.05 0.64 0.16 0.1],'style','popup','string',...
                '*.bmp|*.emf|*.eps|*.jpeg|*.pdf|*.png|*.tiff',...
                'backgroundcolor',[R2 G2 B2],'tag','pripona');
            
            prip=get(findobj('tag','nmin'),'userdata');     
            if isempty(prip)
                set(findobj('tag','pripona'),'value',6);
            else  
                set(findobj('tag','pripona'),'value',prip);
            end
            
            %Text next to the Graphic format file
            g1=uicontrol('Units','normalized','position',[0.3 0.67 0.65 0.13],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Select one of the graphic format files. For a vector output it is recommended to use *.eps graphic file and *.png format for a bitmap output.');
            set(g1,'HorizontalAlignment','left');
            
            %Text Colormap
            uicontrol('Units','normalized','position',[0.05 0.57 0.16 0.05],...
                'style','text','string','Colormap',...
                'backgroundcolor',[R1 G1 B1]);
            
            %Colormap pop-up menu
            uicontrol('Units','normalized','position',[0.05 0.46 0.16 0.1],...
                'style','popup','string',...
                'jet|HSV|hot|cool|spring|summer|autumn|winter|gray|bone|copper|pink|lines',...
                'backgroundcolor',[R2 G2 B2],'tag','colormap');
            set(findobj('tag','colormap'),'value',2);
            
            color=get(findobj('tag','nmax'),'userdata');     
            if isempty(color)
                set(findobj('tag','colormap'),'value',1);
            else  
                set(findobj('tag','colormap'),'value',color);
            end
            
            %Text next to the colormap
            g1=uicontrol('Units','normalized','position',[0.3 0.417 0.65 0.2],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Select a colormap of the output file. Mostly it is recommended to use the jet colormap, which ranges from blue to red, and passes through the colors cyan, yellow, and orange.');
            set(g1,'HorizontalAlignment','left');
            
            %Text Number of colors
            uicontrol('Units','normalized','position',[0.04 0.39 0.19 0.05],...
                'style','text','string','Number of colors',...
                'backgroundcolor',[R1 G1 B1]);
            
            %Value of number of colors
            uicontrol('Units','normalized','position',[0.08 0.32 0.1 0.06],...
                'style','edit','string',15,'backgroundcolor',...
                [R2 G2 B2],'tag','skala');
            
            %Text next to the number of colors
            g1=uicontrol('Units','normalized','position',[0.3 0.31 0.65 0.13],...
                'style','text','backgroundcolor',[R1 G1 B1]);
            set(g1,'HorizontalAlignment','left','string',...
                'Enter a number of colors of the selected colormap. Note that processing time may increase to a several minutes, if a large number of colors has been entered for a large data set.');
            
            ncolor=get(findobj('tag','text_nmin'),'userdata');
            if isempty(ncolor)                
            else
                set(findobj('tag','skala'),'string',ncolor);
            end
            
            %Text DPI 
            uicontrol('Units','normalized','position',[0.08 0.24 0.1 0.05],...
                'style','text','string','DPI','backgroundcolor',[R1 G1 B1]);
            
            %Value of DPI
            uicontrol('Units','normalized','position',[0.08 0.17 0.1 0.06],...
                'style','edit','string',300,'backgroundcolor',[R2 G2 B2],...
                'tag','DPI');
            
            DPI=get(findobj('tag','text_nmax'),'userdata');
            if isempty(ncolor)                
            else
                set(findobj('tag','DPI'),'string',DPI);
            end
            
            %Text next to the DPI
            g1=uicontrol('Units','normalized','position',[0.3 0.03 0.65 0.2],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Enter a value of dots per inch of the output file.');
            set(g1,'HorizontalAlignment','left');
            
            %Panel
            uipanel('Units','normalized','position',[0 0 0.26 1],...
                'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],...
                'fontsize',8,'fontname','cambria');
            
            %OK button
            uicontrol('Units','normalized','position',[0.35 0.02 0.13 0.08],...
                'style','pushbutton','string','OK','Callback','GrafLab OKDDS');
    
            %Close button
            uicontrol('Units','normalized','position',[0.55 0.02 0.13 0.08],...
                'style','pushbutton','string','Close','Callback',...
                'GrafLab CloseDDS');
            
            %setting font to Cambria
            set(get(D,'children'),'fontname','cambria','fontsize',10)
            
        case('fnALFs') %Click on the Computation of fnALFs
            
            %Main window
            F=figure('units','pixels','numbertitle','off','name',...
                'Approach for the computation of fnALFs','color',...
                [R1 G1 B1],'position',[320 150 550 350],'tag',...
                'oknofnALFs','menubar','none'); 
            
            %Radio button group
            u0 = uibuttongroup('visible','on','units','normalized',...
                'Position',[0 0 .26 1],'tag','volbaALFs'); 

            %The first radiobutton
            u1 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.7 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton1'); 
            set(u1,'String','<html>Standard forward<br>column method',...
                'fontname','cambria','fontsize',10);
            note1=uicontrol('Units','normalized','position',...
                [0.3 0.74 0.65 0.2],'style','text','backgroundcolor',...
                [R1 G1 B1]);
            set(note1,'HorizontalAlignment','left','string',...
                'It is recommended to use the standard forward column method for all latitudes up to the maximum degree 1800. However, this method may also be used for the latitudes <0�,56�> and <80�,90�> up to the maximum degree 2190.');
            
            %The second radiobutton
            u2 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.42 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton2'); 
            set(u2,'String',...
                '<html>Modified forward<br>column method<br>combined with<br>Horner큦 scheme',...
                'fontname','cambria','fontsize',10);
            note2=uicontrol('Units','normalized','position',[0.3 0.38 0.65 0.3],...
                'style','text','backgroundcolor',[R1 G1 B1]);
            set(note2,'HorizontalAlignment','left','string',...
                'It is recommended to use the modified forward column method combined with Horner큦 scheme for all latitudes and maximum degrees ranging from 1801 to 2700. This method may also be used for lower degrees than 1801, but cannot be applied for higher degrees than 2700 due to the overflow problem.');            
            
            %The third radiobutton
            u3 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.12 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton3'); 
            set(u3,'String','<html>Extended-range<br>arithmetic',...
                'fontname','cambria','fontsize',10);
            note3=uicontrol('Units','normalized','position',...
                [0.3 0.018 0.65 0.3],'style','text','backgroundcolor',...
                [R1 G1 B1]);
            set(note3,'HorizontalAlignment','left','string',...
                'The extended-range arithmetic approach may be used for all latitudes up to an arbitrary degree essentially.');                        
            
            %OK button
            uicontrol('Units','normalized','position',[0.35 0.02 0.13 0.08],...
                'style','pushbutton','string','OK','Callback','GrafLab OKfnALFs');
    
            %Close button
            uicontrol('Units','normalized','position',[0.55 0.02 0.13 0.08],...
                'style','pushbutton','string','Close','Callback','GrafLab ClosefnALFs');
            
            %The chosen approach for computation of fnALFs
            volbaALFs=get(findobj('tag','fnALFs'),'userdata');
            if isempty(volbaALFs)
                volbaALFs=1;
            end
            
            %Mark the chosen radiobutton
            if volbaALFs==1
                set(u1,'value',1);
            elseif volbaALFs==2
                set(u2,'value',1);
            elseif volbaALFs==3
                set(u3,'value',1);
            end
            
            set(get(F,'children'),'fontname','cambria','fontsize',10);
            
        case('ClosefnALFs') %Click on the Close button in the Computation of fnALFs window
            
            close
            
        case('OKfnALFs') %Click on the OK button in the Computation of fnALFs window
            
            volbaALFs=get(findobj('tag','volbaALFs'),'selectedobject');

            if get(findobj('tag','rbutton1'),'value')==1
                volbaALFs=1;
            elseif get(findobj('tag','rbutton2'),'value')==1
                volbaALFs=2;
            elseif get(findobj('tag','rbutton3'),'value')==1
                volbaALFs=3;
            end
            
            set(findobj('tag','fnALFs'),'userdata',volbaALFs);            
            
            close
            
        case('OKDDS') %Click on the OK button in the Display data settings window
            
            display_data=get(findobj('tag','Display'),'value');
            set(findobj('tag','DDS'),'userdata',display_data);  
            
            prip=get(findobj('tag','pripona'),'value');
            set(findobj('tag','nmin'),'userdata',prip);
            
            color=get(findobj('tag','colormap'),'value');
            set(findobj('tag','nmax'),'userdata',color);

            ncolor=get(findobj('tag','skala'),'string');
            ncolor=str2double(ncolor);
            set(findobj('tag','text_nmin'),'userdata',ncolor);
            
            DPI=get(findobj('tag','DPI'),'string');
            DPI=str2double(DPI);
            set(findobj('tag','text_nmax'),'userdata',DPI);
            
            %Check of the entered value of DPI and number of colors
            if isnan(DPI)==1 || DPI<0
                errordlg('Entered DPI value must be larger than zero.',...
                    'Display data settings');
                error('Entered DPI value must be larger than zero.');
            elseif isnan(ncolor)==1
                errordlg('Entered value of number of colors is not correct.',...
                    'Display data settings');
                error('Entered value of number of colors is not correct.');
            end
            
            if ncolor<2
                errordlg('Value of number of colors must be larger than 1.',...
                    'Display data settings');
                error('Value of number of colors must be larger than 1.');
            end
            
            close
            
        case('CloseDDS') %Click on the Close button in the Display data settings window
            
            close
            
        case('OK') %Click on the OK button in the main GrafLab window 
        
            set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow; 
            
            %Identification of point type selection
            volbagridcheck=get(findobj('tag','gridcheck'),'value'); %Grid
            volbaloadcheck=get(findobj('tag','loaddata'),'value'); %Load
            volbadiskcheck=get(findobj('tag','diskcheck'),'value'); %Diskretne body
            
            %Error messages for point type selection
            if volbagridcheck==0 && volbadiskcheck==0 && volbaloadcheck==0
                errordlg('Please select a point type.','Point type selection error');
                error('Please select a point type.');
            elseif volbagridcheck==1 && volbadiskcheck==1 && volbaloadcheck==1
                errordlg('Please select one point type only.',...
                    'Error in point type selection');
                error('Please select one point type only.');
            elseif volbagridcheck==1 && volbadiskcheck==1
                errordlg('Please select one point type only.',...
                    'Error in point type selection');
                error('Please select one point type only.');
            elseif volbadiskcheck==1 && volbaloadcheck==1
                errordlg('Please select one point type only.',...
                    'Error in point type selection');
                error('Please select one point type only.');
            elseif volbagridcheck==1 && volbaloadcheck==1
                errordlg('Please select one point type only.',...
                    'Error in point type selection');
                error('Please select one point type only.');
            end            
                            
            %Ellipsoid
            if get(findobj('tag','ell'),'value') == 1; %GRS80              
                GMEl=3986005*10^8; %Geocentric gravitational constant of GRS80
                aEl=6378137; %Semimajor axis of GRS80
                eEl=sqrt(0.00669438002290342); %First eccentricity of GRS80
                omegaEl=7292115*10^-11; %Angular velocity of GRS80
                CEl_20=-108263*10^-8/sqrt(5); %Fully normalized C_20 of GRS80
            elseif get(findobj('tag','ell'),'value') == 2; %WGS84
                GMEl=3986004.418*10^8; %Geocentric gravitational constant of WGS84
                aEl=6378137; %Semimajor axis of WGS84
                fEl=1/298.257223563; %Flattening of WGS84
                omegaEl=7292115*10^-11; %Angular velocity of WGS84
                CEl_20=-0.484166774985*10^-3; %Fully normalized C_20 of WGS84
                eEl=sqrt(fEl*(2-fEl)); %First eccentricity of WGS84
            end         
            
            %GM a R of GGM
            GM=str2double(get(findobj('tag','GM'),'string'));            
            R=str2double(get(findobj('tag','R'),'string'));
            if GM<=0
                errordlg('Value of GM must be larger than zero.',...
                    'Error in geopotential model and reference system selection')
                error('Value of GM must be larger than zero.')
            elseif R<=0
                errordlg('R must be larger than zero.',...
                    'Error in geopotential model and reference system selection')
                error('R value must be larger than zero.')
            end            
                                           
            %Selection of functionals of the geopotential
            volbapar1=get(findobj('tag','P1'),'value'); %ID of the first functional
            volbapar2=get(findobj('tag','P2'),'value'); %ID of the second functional
            volbapar3=get(findobj('tag','P3'),'value'); %ID of the third functional
            volbapar4=get(findobj('tag','P4'),'value'); %ID of the fourth functional
            volbapar=[volbapar1;volbapar2;volbapar3;volbapar4];
            pocetpar=length(find(volbapar>1));

            %Error messages for selection of functionals of the geopotential
            if pocetpar==1
                if volbapar1>1 && volbapar2==1 && volbapar3==1 && volbapar4==1
                else
                    errordlg('Please select the functional of the geopotential in the first pop-up menu.',...
                        'Calculated parameters and output selection');
                    error('Please select the functional of the geopotential in the first pop-up menu.');
                end
            elseif pocetpar==2
                if volbapar1>1 && volbapar2>1 && volbapar3==1 && volbapar4==1
                    if volbapar1 == volbapar2
                        errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                            'Calculated parameters and output selection');
                        error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                    end
                else
                    errordlg('Please select the functionals of the geopotential in the first and the second pop-up menu.',...
                        'Calculated parameters and output selection');
                    error('Please select the functionals of the geopotential in the first and the second pop-up menu.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 %If geoid or height anomaly is to be computed
                    errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                        'Calculated parameters and output selection');
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            elseif pocetpar==3
                if volbapar1>1 && volbapar2>1 && volbapar3>1 && volbapar4==1
                    if (volbapar1 == volbapar2) || (volbapar1 == volbapar3) || (volbapar2 == volbapar3)
                        errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                            'Calculated parameters and output selection');
                        error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                    end
                else
                    errordlg('Please select the functionals of the geopotential in the first, second and third pop-up menu.',...
                        'Calculated parameters and output selection');
                    error('Please select the functionals of the geopotential in the first, second and third pop-up menu.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 || length(nonzeros(volbapar==10 | volbapar==23))==2 %If geoid or height anomaly is to be computed
                    errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                        'Calculated parameters and output selection');
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            elseif pocetpar==4
                if (volbapar1 == volbapar2) || (volbapar1 == volbapar3) || (volbapar1 == volbapar4) || (volbapar2 == volbapar3) || (volbapar2 == volbapar4) || (volbapar3 == volbapar4)
                    errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                        'Calculated parameters and output selection');
                    error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 || length(nonzeros(volbapar==10 | volbapar==23))==2 || length(nonzeros(volbapar==10 | volbapar==23))==3 %If geoid or height anomaly is to be computed
                    errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                        'Calculated parameters and output selection');
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            end
            
            %Identification of type of the input coordinates
            coord=get(findobj('tag','rbutton2coord'),'value');

            if volbagridcheck==1
                if coord==1 %Entered spherical coordinates                                                
                    if any(volbapar==10) || any(volbapar==23) 
                        errordlg('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. Grid that has been entered in the spherical coordinates (constant value of the radius r) does not refer to this surface. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates a set value in the array Height above the reference surface (m) to zero.','Calculated parameters and output selection')
                        error('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. Grid that has been entered in the spherical coordinates (constant value of the radius r) does not refer to this surface. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates a set value in the array Height above the reference surface (m) to zero.')
                    end
                end
            elseif volbadiskcheck==1
                if coord==1 %Entered spherical coordinates                                                
                    if any(volbapar==10) || any(volbapar==23) 
                        errordlg('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates a set value(s) in the array Ellipsoidal height/Spherical radius (m) to zero.','Calculated parameters and output selection')
                        error('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates a set value(s) or in the array Ellipsoidal height/Spherical radius (m) to zero.')
                    end
                end
            end
            
            if volbapar1==1 && volbapar2==1 && volbapar3==1 && volbapar4==1
                errordlg('Please choose at least one functional of the geopotential.',...
                    'Error in calculated paramters and output selection')
                error('Please choose at least one functional of the geopotential.')   
            end 
            
            %Identification of computation of fnALFs approach
            volbaALFs=get(findobj('tag','fnALFs'),'userdata');
            if isempty(volbaALFs)
                volbaALFs=1;
            end
            
            %Identification of Display data
            display_data=get(findobj('tag','DDS'),'userdata');
            if isempty(display_data)
                display_data=0;
            end
            
            %Identification of computation of commission error
            STD=get(findobj('tag','STD'),'value');
            
            %Identification of output folderu and file
            outname=get(findobj('tag','R_text'),'userdata');
            outadresar=get(findobj('tag','ell_text'),'userdata');
            if isempty(outname)
                set(findobj('tag','hlasky'),'string','Select output folder and output file',...
                       'fontsize',8,'foregroundcolor','k'); drawnow;
                
                warn2=warndlg('Output folder and output file were not specified. Click OK and then select an output folder and output file.');
                waitfor(warn2);
                
                [outname,outadresar]=uiputfile('*.*');
                if outname==0
                    errordlg('Output folder and output file must be specified!');
                    
                    set(findobj('tag','hlasky'),'string','',...
                       'fontsize',8,'foregroundcolor','k'); drawnow;
                   
                    error('Output folder and output file must be specified!');
                else
                    if find(outname=='.')>0
                        outname=outname(1:(find(outname=='.')-1));
                    end
                end
                
                set(findobj('tag','hlasky'),'string','',...
                       'fontsize',8,'foregroundcolor','k'); drawnow;
            end
            
            %Check of the entered values of DPI a number of colors
            if display_data==1
                ncolor=get(findobj('tag','text_nmin'),'userdata');
                DPI=get(findobj('tag','text_nmax'),'userdata');
                if isnan(DPI)==1 || DPI<0
                    errordlg('Entered DPI value must be larger than zero.',...
                        'Display data settings');
                    error('Entered DPI value must be larger than zero.');
                elseif isnan(ncolor)==1 || ncolor<2
                    errordlg('Value of number of colors must be larger than 1.',....
                        'Display data settings');
                    error('Value of number of colors must be larger than 1.');
                end
            end
            
            %Check, if computation of tensors in the LNOF has been selected (not allowed)
            if any(volbapar==8) || any(volbapar==9) || any(volbapar==14) || any(volbapar==15)
                if volbaALFs==2 && volbagridcheck==1 && STD==0
                    errordlg('The following functionals of the geopotential: Disturbing_tensor_Txx_Tyy_Tzz, Disturbing_tensor_Txy_Txz_Tyz, Gravitational_tensor_Vxx_Vyy_Vzz and Gravitational_tensor_Vxy_Vxz_Vyz cannot be computed on a grid using the modified forward column method combined with Horner큦 scheme. In this case, for a high maximum degree (1800 and higher), it is recommended to use the extended-range arithmetic approach. However, for the point-wise computation of the above mentioned functionals, the modified forward column method combined with Horner큦 scheme may be applied.',...
                        'Calculated parameters and output selection');
                    error('The following functionals of the geopotential: Disturbing_tensor_Txx_Tyy_Tzz, Disturbing_tensor_Txy_Txz_Tyz, Gravitational_tensor_Vxx_Vyy_Vzz and Gravitational_tensor_Vxy_Vxz_Vyz cannot be computed on a grid using the modified forward column method combined with Horner큦 scheme. In this case, for a high maximum degree (1800 and higher), it is recommended to use the extended-range arithmetic approach. However, for the point-wise computation of the above mentioned functionals, the modified forward column method combined with Horner큦 scheme may be applied.');
                end
            end
             
            %Check, if display data in point-wise or load data approach has been
            %selected (not allowed)
            if display_data==1 && volbagridcheck~=1
                warn4=warndlg('Only data computed on a grid may be displayed. After clicking OK, the computation will start, although the computed data will not be displayed.');
                waitfor(warn4);
                display_data=0;
            end           
            
            if STD==0
              
                tic %Start clock to measure computation time
                
                %Loading DMR
                if get(findobj('tag','P1'),'value')==10 || get(findobj('tag','P2'),'value')==10 || get(findobj('tag','P3'),'value')==10 || get(findobj('tag','P4'),'value')==10 || get(findobj('tag','P1'),'value')==23 || get(findobj('tag','P2'),'value')==23 || get(findobj('tag','P3'),'value')==23 || get(findobj('tag','P4'),'value')==23
                   set(findobj('tag','hlasky'),'string','Please select DTM file',...
                       'fontsize',8,'foregroundcolor','k'); drawnow;
                   [loadnameDMR,loadadresarDMR]=uigetfile('*.*');

                   if loadnameDMR==0
                       errordlg('In order to compute Geoid_undulation/Height_anomaly, DTM file must be imported!',...
                           'Error in geopotential model and reference system selection');
                       error('In order to compute Geoid_undulation/Height_anomaly, DTM file must be imported!')
                   else
                       set(findobj('tag','hlasky'),'string',...
                           'Loading DTM file...','fontsize',8,...
                           'foregroundcolor','k'); drawnow;

                       DMR=importdata([loadadresarDMR,loadnameDMR]);
                       HC=DMR(:,3);
                       HS=DMR(:,4);

                       clear DMR
                       
                       set(findobj('tag','hlasky'),'string','',...
                           'foregroundcolor','k'); drawnow;
                   end
                end                                      

                %Loading GGM                              
                GGMname=get(findobj('tag','R'),'userdata');
                GGMadresar=get(findobj('tag','ell'),'userdata');

                if isempty(GGMname) %Error message, if GGM file has not been imported
                    errordlg('Please input geopotential model file.',...
                        'Error in point type selection');
                    error('Please input geopotential model file.')
                end

                set(findobj('tag','hlasky'),'string',...
                        'Loading GGM file...','fontsize',8,...
                        'foregroundcolor','k'); drawnow;

                %Loading GGM
                if strcmp(GGMname(end-3:end),'.gfc') %Input data in ICGEM format 
                    fGGMid=fopen([GGMadresar,GGMname]);

                    while(~feof(fGGMid))                       
                        s=fgetl(fGGMid); 
                        if strncmpi(s,'end_of_head',11)
                            break
                        end
                    end

                    GGM=textscan(fGGMid,'%s%f%f%f%f%f%f');

                    fclose(fGGMid);

                    GGM=GGM(2:end);
                    GGM=cell2mat(GGM);
                else %Input data in txt file or mat file
                    GGM=importdata([GGMadresar,GGMname]);
                end
                
                stupen=GGM(:,1);
                rad=GGM(:,2);
                C=GGM(:,3);
                S=GGM(:,4);

                del00=find(stupen==0 & rad==0);
                C(del00)=[];
                S(del00)=[];
                stupen(del00)=[];
                rad(del00)=[];

                del10=find(stupen==1 & rad==0);
                C(del10)=[];
                S(del10)=[];
                stupen(del10)=[];
                rad(del10)=[];

                del11=find(stupen==1 & rad==1);
                C(del11)=[];
                S(del11)=[];
                stupen(del11)=[];
                rad(del11)=[];

                %Identification of GGM file format
                if stupen(1)==2 && stupen(2)==2 && rad(1)==0 && rad(2)==1
                    radenie=0;
                elseif stupen(1)==2 && stupen(2)==3 && rad(1)==0 && rad(2)==0
                    radenie=1;                   
                else
                    errordlg('Wrong format of the inputted GGM file',...
                        'Geopotential model and reference system selection');
                    error('Wrong format of the inputted GGM file')
                end                

                set(findobj('tag','hlasky'),'string',...
                        '','fontsize',8,'foregroundcolor','k'); drawnow;
                
                %Value of nmin and error messages
                nmin=str2double(get(findobj('tag','nmin'),'string'));
                nmaxGGM=max(GGM(:,1));
                if nmin<0 
                    errordlg('Value of nmin cannot be negative.',...
                        'Error in geopotential model and reference system selection')
                    error('Value of nmin cannot be negative.')
                elseif nmin>nmaxGGM
                    errordlg('Value of nmin exceedes nmax value of GGM.',...
                        'Error in geopotential model and reference system selection')
                    error('Value of nmin exceedes nmax value of GGM.')
                end
                if isnan(nmin)==1
                        errordlg('Please input the nmin value.',...
                            'Error in geopotential model and reference system selection')
                        error('Please input the nmin value.')
                end
                if rem(nmin,1)~=0
                    errordlg('Value of nmin must be an integer.',...
                        'Error in geopotential model and reference system selection')
                    error('Value of nmin must be an integer.')
                end
                
                %Value of nmax and error messages   
                if get(findobj('tag','use'),'value')==1
                    nmax=nmaxGGM;
                else
                    nmax=str2double(get(findobj('tag','nmax'),'string'));

                    if nmax>nmaxGGM
                        errordlg('Entered value of nmax exceedes nmax value of GGM.',...
                            'Error in geopotential model and reference system selection')
                        error('Entered value of nmax exceedes nmax value of GGM.')
                    elseif nmax<2
                        errordlg('Value of nmax must be at least 2.',...
                            'Error in geopotential model and reference system selection')
                        error('Value of nmax must be at least 2.')
                    elseif nmin>nmax
                        errordlg('Value of nmin cannot be larger than nmax value.',...
                            'Error in geopotential model and reference system selection')
                        error('Value of nmin cannot be larger than nmax value.')
                    end
                    if isnan(nmax)==1
                        errordlg('Please input the nmax value.',...
                            'Error in geopotential model and reference system selection')
                        error('Please input the nmax value.')
                    end
                    
                    if rem(nmax,1)~=0
                        errordlg('Value of nmax must be an integer.',...
                            'Error in geopotential model and reference system selection')
                        error('Value of nmax must be an integer.')
                    end
                end  
                
                if nmin==0 %Zero degree term for gravitational and gravity funtionals
                    nmin=2;
                    nmin0=1; %Logical 1
                    nmin1=0; %Logical 0
                else
                    nmin0=0; %Logical 0
                    nmin1=0; %Logical 0
                end
                
                if nmin==1
                    nmin0=0; %Logical 0
                    nmin1=1; %Logical 1
                    nmin=2;                 
                end
                
                if nmin>2
                    if any(volbapar==10) || any(volbapar==23)
                        errordlg('Geoid_undulation and Height_anomaly cannot be computed if nmin>2. Please use the functional: Height_anomaly_ell and compute these functionals in the ellipsoidal coordinates with the array Height above the reference surface (m) (or Ellipsoidal height/Spherical radius (m)) set to zero, or compute Geoid_undulation/Height_anomaly twice. The first time with nmax equal to your value of nmax and the second time with nmax=nmin-1. In both cases set nmin to zero. Differece between these to values is your required Geoid_undulation/Height_anomaly.','Calculated parameters and output selection');
                        error('Geoid_undulation and Height_anomaly cannot be computed if nmin>2. Please use the functional: Height_anomaly_ell and compute these functionals in the ellipsoidal coordinates with the array Height above the reference surface (m) (or Ellipsoidal height/Spherical radius (m)) set to zero, or compute Geoid_undulation/Height_anomaly twice. The first time with nmax equal to your value of nmax and the second time with nmax=nmin-1. In both cases set nmin to zero. Differece between these to values is your required Geoid_undulation/Height_anomaly.');
                    end
                    
                    if any(volbapar==9) || any(volbapar==15)
                        errordlg('The following functionals of the geopotential: Disturbing_tensor_Txy_Txz_Tyz and Gravitational_tensor_Vxy_Vxz_Vyz cannot be computed if nmin>2.',...
                            'Calculated parameters and output selection');
                        error('The following functionals of the geopotential: Disturbing_tensor_Txy_Txz_Tyz and Gravitational_tensor_Vxy_Vxz_Vyz cannot be computed if nmin>2.');
                    end
                    
                    if any(volbapar==8) || any(volbapar==14) 
                        LNOFnmin=1; %Logical 1
                        
                        if any(volbapar~=8 & volbapar~=9 & volbapar~=14 & volbapar~=15 & volbapar~=1)                                              
                            errordlg('The following functionals of the geopotential: Disturbing_tensor_Txx_Tyy_Tzz and Gravitational_tensor_Vxx_Vyy_Vzz cannot be computed simultaneously with other functionals if nmin>2.',...
                            'Calculated parameters and output selection');
                            error('The following functionals of the geopotential: Disturbing_tensor_Txx_Tyy_Tzz and Gravitational_tensor_Vxx_Vyy_Vzz cannot be computed simultaneously with other functionals if nmin>2.');
                        end
                    else
                        LNOFnmin=0; %Logical 0
                    end
                end                              
                
                %% Computation of functionals of the geopotential on a regular grid           

                if volbagridcheck==1

                    %Entered coordinates of the grid
                    fimin=str2double(get(findobj('tag','fimin'),'string')); 
                    fistep=str2double(get(findobj('tag','fistep'),'string'));
                    fimax=str2double(get(findobj('tag','fimax'),'string'));
                    lambdamin=str2double(get(findobj('tag','lambdamin'),'string')); 
                    lambdastep=str2double(get(findobj('tag','lambdastep'),'string'));
                    lambdamax=str2double(get(findobj('tag','lambdamax'),'string'));
                    h=str2double(get(findobj('tag','hgrid'),'string'));

                    %Check of the entered coordinates
                    if isnan(fimin) || isnan(fistep) || isnan(fimax) || isnan(lambdamin) || isnan(lambdastep) || isnan(lambdamax) || isnan(h)
                        errordlg('Entered grid is not correct.',...
                            'Error in point type selection');
                        error('Entered grid is not correct.');       
                    end

                    if fimin>fimax
                        errordlg('Value of Lat. min must be smaller than the Lat. max value.',...
                            'Error in point type selection');
                        error('Value of Lat. min must be smaller than the Lat. max value.'); 
                    elseif fistep<=0
                        errordlg('Value of Lat. step must be larger than zero.',...
                            'Error in point type selection');
                        error('Value of Lat. step must be larger than zero.'); 
                    elseif lambdamin>lambdamax
                        errordlg('Value of Lon. min must be smaller than Lon. max value.',...
                            'Error in point type selection');
                        error('Value of Lon. min must be smaller than Lon. max value.');
                    elseif lambdastep<=0
                        errordlg('Lon. step must be larger than zero.',...
                            'Error in point type selection');
                        error('Lon. step must be larger than zero.'); 
                    end

                    if fimin>90 || fimin<-90
                        errordlg('Value of Lat. min must be within the interval <-90�,90�>.',...
                            'Error in point type selection');
                        error('Value of Lat. min must be within the interval <-90�,90�>.');
                    end
                    if fimax>90 || fimax<-90
                        errordlg('Value of Lat. max must be within the interval <-90�,90�>.',...
                            'Error in point type selection');
                        error('Value of Lat. max must be within the interval <-90�,90�>.');
                    end
                    if lambdamin>360 || lambdamin<-180
                        errordlg('Value of Lon. min must be within the interval <-180�,180�> or <0�,360�>.',...
                            'Error in point type selection');
                        error('Value of Lon. min must be within the interval <-180�,180�> or <0�,360�>.');
                    end
                    if lambdamax>360 || lambdamax<-180
                        errordlg('Value of Lon. max must be within the interval <-180�,180�> or <0�,360�>.',...
                            'Error in point type selection');
                        error('Value of Lon. max must be within the interval <-180�,180�> or <0�,360�>.');
                    end
                    if (lambdamax-lambdamin)>360
                        errordlg('Longitude must be in the range <-180�,180�> or <0�,360�>.',...
                            'Error in point type selection');
                        error('Longitude must be in the range <-180�,180�> or <0�,360�>.');
                    end
                    
                    %Vectors phi and lambda in one longitude and latitude
                    %parallel, respectively
                    fi=(fimin:fistep:fimax)';
                    lambda=(lambdamin:lambdastep:lambdamax)';  
                    fi=deg2rad(fi(:));
                    lambda=deg2rad(lambda(:));

                    %Grid, which is to be displayed has to have at least two
                    %points in latitude parallels and at least two points in
                    %longitude parallels.
                    if display_data==1
                        if length(fi)<2 || length(lambda)<2
                            warn3=warndlg('In order to display the computed data on the grid, it has to contains at least two different points in one latitude parallel and two different points in one longitude longitude. After clicking OK, the computation will start, although the data will not be displayed.');
                            waitfor(warn3);
                            display_data=0;
                        end
                    end
                    
                    if coord==1 %Entered spherical coordinates                       
                        %Spherical radius
                        r=(R+h)*ones(length(fi),1);
                        hsph=h;
                        
                        %Spherical latitude
                        fiG=fi;  
                        
                        %Transform spherical latitude into the ellipsoidal
                        %latitude
                        [X Y Z]=sph2cart(0*fiG,fiG,r);
                        [fi , ~, h]=ecef2geodetic(X,Y,Z,[aEl eEl]');

                        clear X Y Z
                    elseif coord==0 %Entered ellipsoidal coordinates
                        %Trasformation of (fi, lambda, h) into (X, Y, Z)
                        [X,Y,Z]=geodetic2ecef(fi,0*zeros(length(fi),1),h*ones(length(fi),1),[aEl eEl]');  
                        r=sqrt(X.*X+Y.*Y+Z.*Z); %Radius

                        %Spherical latitude
                        fiG=atan(Z./sqrt(X.*X+Y.*Y)); 
                        
                        clear X Y Z 
                    end

                    %Computation of the coefficients C2,0; C4,0; C6,0;
                    %C8,0; C10,0 of the selected ellipsoid
                    CEl=zeros(length(C),1);
                    for n=1:5
                        CEl(2*n==stupen & rad==0,1)=((-1)^n*(3*eEl^(2*n))/((2*n+1)*(2*n+3)*sqrt(4*n+1))*(1-n-5^(3/2)*n*CEl_20/eEl^2)).*(aEl./R).^(2*n).*(GMEl/GM);
                    end                                

                    if any(volbapar==11) || any(volbapar==12) || any(volbapar==13) || any(volbapar==14) || any(volbapar==15) || any(volbapar==16) || any(volbapar==17) || any(volbapar==18) || any(volbapar==20) || any(volbapar==25)
                        grav=1;
                    else
                        grav=0;
                    end

                    if  any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==5) || any(volbapar==6) || any(volbapar==7) || any(volbapar==8) || any(volbapar==9) || any(volbapar==10) || any(volbapar==19) || any(volbapar==21) || any(volbapar==22) || any(volbapar==23) || any(volbapar==24)
                        por=1;
                        deltaC=C-CEl;
                    else
                        por=0;
                    end
 
                    if any(volbapar==20)
                        normal=1;
                    else
                        normal=0;
                    end

                    clear GGM stupen rad                   
                    if normal==0
                        clear CEl
                    end

                    %Initialization
                    eta=0; ksi=0; Theta=0; T=0; T_rr=0; Trr=0; Trf=0; Trl=0; Tff=0;
                    Tfl=0; Tll=0; Tzz=0; Txx=0; Tyy=0; N=0; V=0; Vrr=0; Vrf=0; Vrl=0; Vff=0;
                    Vfl=0; Vll=0; g=0; g_sa=0; W=0; anomalia_sa=0; porucha=0; 
                    porucha_sa=0; zetaEl=0; zeta=0; Wrr=0; Wr=0; Wfi=0; Wlambda=0;
                    Ur=0; Ufi=0; N1c=0; N2c=0; H=0;

                    %If geoid/quasigeoid is to be computed
                    if any(volbapar==10) || any(volbapar==23)
                        geoid=1;
                        if h~=0
                            errordlg('In order to compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                                'Error in point type selection');
                            error('In order to compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                        end
                    else
                        geoid=0;
                    end

                    %Indices of the spherical harmonic coefficients
                    if radenie==0
                        index=zeros(nmax-1,1);
                        index(1,1)=1;
                        for i=1:(nmax-2)
                            index(i+1,1)=index(i,1)+2+i;
                        end
                    elseif radenie==1
                        z=((nmaxGGM+1)*(nmaxGGM+2)-6)/2-(sum(nmaxGGM-(nmaxGGM:-1:(nmax-1)))-1);
                        k=z;
                    end

                    %Initialization of the matrices and vectors for the computation of fnALFs
                    Pnm=zeros(length(fiG),nmax+1);
                    q=(R./r);
                    q2=(R./r).^2;
                    u=cos(fiG);
                    t=sin(fiG);
                    
                    %Initialization for extended-range arithmetic approach
                    if volbaALFs==3
                                               
                        bit=mexext; %Bit version of Matlab
                        bit=bit(end-1:end);
                        bit=str2double(bit);

                        nmax23=nmax*2+3;
                        rr=zeros(nmax23,1); ri=rr;
                        dd=zeros(nmax,1); am=dd; bm=am;

                        m1=1:nmax23;
                        rr(m1)=sqrt(m1);
                        ri(m1)=1./rr;
                        m2=1:nmax;
                        dd(m2)=rr(2*m2+3).*ri(2*m2+2);

                        IND=960;
                        BIG=2^IND;
                        BIGI=2^(-IND);
                        BIGS=2^(IND/2);
                        BIGSI=2^(-IND/2);
                        ROOT3=1.732050807568877;
                        
                        if bit==32
                            pm=am;
                            ps1=zeros(length(fiG),nmax); 
                            ips1=ps1;
                            x=ROOT3*u.*q;
                            ix=zeros(size(x));
                            ps1(:,1)=x;
                            ips1(:,1)=ix;
                            for m3=2:nmax
                                x=(dd(m3-1)*u).*x.*q;
                                y=abs(x);
                                iy=y>=BIGS;
                                if any(iy)
                                    x(iy)=x(iy)*BIGI;
                                    ix(iy)=ix(iy)+1;
                                end
                                iy=y<BIGSI;
                                if any(iy)
                                    x(iy)=x(iy)*BIG;
                                    ix(iy)=ix(iy)-1;
                                end
                                ps1(:,m3)=x;
                                ips1(:,m3)=ix;
                            end
                        elseif bit==64
                            zzb=zeros(length(fiG),1);
                            ps1b=zeros(length(fiG),nmax); 
                            ips1b=ps1b;
                            xb=ROOT3*u.*q;
                            ixb=zeros(size(xb));
                            ps1b(:,1)=xb;
                            ips1b(:,1)=ixb;
                            for m3=2:nmax
                                xb=(dd(m3-1)*u).*xb.*q;
                                yb=abs(xb);
                                iyb=yb>=BIGS;
                                if any(iyb)
                                    xb(iyb)=xb(iyb)*BIGI;
                                    ixb(iyb)=ixb(iyb)+1;
                                end
                                iyb=yb<BIGSI;
                                if any(iyb)
                                    xb(iyb)=xb(iyb)*BIG;
                                    ixb(iyb)=ixb(iyb)-1;
                                end
                                ps1b(:,m3)=xb;
                                ips1b(:,m3)=ixb;
                            end
                        end
                        
                        clear dd
                    end

                    %Initialization of the matrices and vectors for the 
                    %computation of the first-order derivatives of fnALFs
                    if any(volbapar==3) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16)  || any(volbapar==20)
                        dALFs=1;
                        dPnm=zeros(length(fi),nmax+1);
                        qu=q./u;
                        tu=t./u;
                        
                        %Treatment of the dPnm singularity
                        singdPnm=fi==pi/2 | fi==-pi/2;
                    else
                        dALFs=0;
                    end   
                    
                    %Initialization of the matrices and vectors for the 
                    %computation of the second-order derivatives of fnALFs
                    if any(volbapar==6) || any(volbapar==12)
                        ddALFs=1;
                        ddPnm=zeros(length(fi),nmax+1);
                        
                        %Treatment of the ddPnm singularity
                        singddPnm=fi==pi/2 | fi==-pi/2;
                    else
                        ddALFs=0;
                    end   
                    
                    %Status line
                    progressbar=findobj('tag','hlasky');

                    %% Summation over m
                    for m=nmax:-1:0

                        %Update of the progress bar
                        if rem(m,10)==0
                            set(progressbar,'string',...
                                sprintf('Progress: m = %5.0d',m),...
                                'fontsize',8); drawnow;
                        end

                        %Selection of the spherical harmonic coefficients of order m
                        %======================================================
                        if m<2  
                            if radenie==0                          
                                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                                   Cm=C(index+m);
                                end

                                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                                    deltaCm=deltaC(index+m);
                                end

                                if normal==1 %C's spherical harmonic coefficients for the functionals with the normal field
                                    if m==0
                                        CElm=CEl(index+m);
                                    end
                                end

                                if geoid==1
                                    if m==1
                                        HCm=HC([3;index+m+3]); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                        HSm=HS([3;index+m+3]); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                    elseif m==0
                                        HCm=HC([1;2;index+m+3]); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                        HSm=HS([1;2;index+m+3]); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                    end
                                end
                                
                                Sm=S(index+m);
                            elseif radenie==1
                                z=z+1;

                                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                                    Cm=C(z:k);
                                end

                                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                                    deltaCm=deltaC(z:k);
                                end

                                if normal==1 %C's spherical harmonic coefficients for the functionals with the normal field
                                    if m==0
                                        CElm=CEl(z:k);
                                    end
                                end
                                
                                if geoid==1
                                    if m==1
                                        HCm=HC((z-1):k);
                                        HSm=HS((z-1):k);
                                    elseif m==0
                                        HCm=HC((z-1):k);
                                        HSm=HS((z-1):k);
                                    end
                                end
                                
                                Sm=S(z:k);
                                
                                zT=z;
                                kT=k;
                                
                                z=z-nmaxGGM+m-1;
                                k=z+nmax-m; 
                            end
                        else
                            if radenie==0
                                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                                    Cm=C(index((m-1):end)+m);
                                end

                                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                                    deltaCm=deltaC(index((m-1):end)+m);
                                end
                                
                                if geoid==1
                                    HCm=HC(index((m-1):end)+m+3); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                    HSm=HS(index((m-1):end)+m+3); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                end

                                Sm=S(index((m-1):end)+m);
                            elseif radenie==1                              
                                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                                    Cm=C(z:k);
                                end

                                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                                    deltaCm=deltaC(z:k);
                                end
                                
                                if geoid==1
                                    HCm=HC(z:k);
                                    HSm=HS(z:k);
                                end

                                Sm=S(z:k);
                                
                                zT=z;
                                kT=k;
                                                               
                                z=z-nmaxGGM+m-2;
                                k=z+nmax-m+1;                                                                                                             
                            end
                        end
                        %======================================================

                        %% Computation of modified fnALFs
                        if volbaALFs==1 %Standard forward column method
                            if m==0
                                Pnm(:,1)=1;
                            elseif m==1                    
                                Pnm(:,1)=sqrt(3)*u.*q;  
                            elseif m>1                            
                                i=2*(2:m);
                                i1=sqrt((i+ones(size(i)))./i);
                                Pnm(:,1)=u.^m*sqrt(3)*prod(i1).*q.^m;
                            end

                            if m==nmax
                            elseif m<=(nmax-1)
                                n=m+1;
                                anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                            end

                            if m<(nmax-1)
                                j=3;
                                for n=m+2:nmax                            
                                    anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                    bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                                    Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                                    j=j+1;
                                end
                            end

                        elseif volbaALFs==2 %Modified forward column method
                            if m==0
                                Pnm(:,1)=1e-280;
                            elseif m==1

                                Pnm(:,1)=sqrt(3)*q*1e-280;  
                            elseif m>1                            
                                i=2*(2:m);
                                i1=sqrt((i+ones(size(i)))./i);
                                Pnm(:,1)=sqrt(3)*prod(i1)*(q.^m)*1e-280;
                            end

                            if m==nmax
                            elseif m<=(nmax-1)
                                n=m+1;
                                anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                            end

                            if m<(nmax-1)
                                j=3;
                                for n=m+2:nmax                            
                                    anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                    bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                                    Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                                    j=j+1;
                                end
                            end

                        elseif volbaALFs==3 %Extended-range arithmetic 
                            if bit==32 %32 bit version of Matlab
                     
                                am(m+1)=rr(2*m+3);
                                for n=m+2:nmax
                                    ww=rr(2*n+1)*ri(n-m)*ri(n+m);
                                    am(n)=rr(2*n-1)*ww;
                                    bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*ww;
                                end

                                if m~=0
                                    for i=1:length(fiG) 
                                        x=ps1(i,m);
                                        ix=ips1(i,m);

                                        if ix==0
                                            pm(m)=x;
                                        elseif ix<-1
                                            pm(m)=0;  
                                        elseif ix<0
                                            pm(m)=x*BIGI;
                                        else
                                            pm(m)=x*BIG;
                                        end

                                        if m==nmax
                                            Pnm(i,1:(nmax-m+1))=pm(m:end);
                                            continue;
                                        end

                                        y=x;
                                        iy=ix;
                                        x=(am(m+1)*t(i)*q(i))*y;
                                        ix=iy;
                                        w=abs(x);

                                        if w>=BIGS
                                            x=x*BIGI;
                                            ix=ix+1;
                                        elseif w<BIGSI
                                            x=x*BIG;
                                            ix=ix-1;
                                        end

                                        if ix==0
                                            pm(m+1)=x;
                                        elseif ix<-1  
                                            pm(m+1)=0;    
                                        elseif ix<0
                                            pm(m+1)=x*BIGI;
                                        else
                                            pm(m+1)=x*BIG;
                                        end

                                        for n=m+2:nmax 
                                            id=ix-iy;

                                            if id==0
                                                zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*y;
                                                iz=ix;
                                            elseif id==1
                                                zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*(y*BIGI);
                                                iz=ix;
                                            elseif id==-1
                                                zz=(am(n)*t(i)*q(i))*(x*BIGI)-bm(n)*q2(i)*y;
                                                iz=iy;
                                            elseif id>1
                                                zz=(am(n)*t(i)*q(i))*x;
                                                iz=ix;
                                            else
                                                zz=-bm(n)*q2(i)*y;
                                                iz=iy;
                                            end

                                            w=abs(zz);

                                            if w>=BIGS
                                                zz=zz*BIGI;
                                                iz=iz+1;
                                            elseif w<BIGSI
                                                zz=zz*BIG;
                                                iz=iz-1;
                                            end

                                            if iz==0
                                                pm(n)=zz;
                                            elseif iz<-1
                                                pm(n)=0;     
                                            elseif iz<0
                                                pm(n)=zz*BIGI;
                                            else
                                                pm(n)=zz*BIG;
                                            end

                                            y=x;
                                            iy=ix;
                                            x=zz;
                                            ix=iz;                                           
                                        end

                                        Pnm(i,1:(nmax-m+1))=pm(m:end);  
                                    end

                                elseif m==0
                                    Pnm(:,1)=1;
                                    Pnm(:,2)=sqrt(3)*t.*q;

                                    for i=2:nmax
                                        Pnm(:,i+1)=Pnm(:,i).*sqrt((2*i+1)*(2*i-1))./i.*t.*q-q2.*Pnm(:,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                                    end

                                    clear rr ri am bm pm ps1 ips1 m1 m2 ...
                                        dd ix x y iy w iz zz
                                end
                            elseif bit==64 %64 bit version of Matlab
                                
                                am(m+1)=rr(2*m+3);
                                for n=m+2:nmax
                                    ww=rr(2*n+1)*ri(n-m)*ri(n+m);
                                    am(n)=rr(2*n-1)*ww;
                                    bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*ww;
                                end
                                
                                if m==0 %Zonal modified fnALFs
                                    Pnm(:,1)=1;
                                    Pnm(:,2)=sqrt(3)*t.*q;

                                    for i=2:nmax
                                        Pnm(:,i+1)=Pnm(:,i).*sqrt((2*i+1)*(2*i-1))./i.*t.*q-q2.*Pnm(:,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                                    end

                                    clear rr ri am bm ps1b ips1b m1 m2 dd ...
                                        ixb xb yb iyb wb izb zzb pmxb pm0b pmxBIGIb ...
                                        pmxBIGb wb wBIGSb wBIGSIb pm1xb pm10b ...
                                        pm1xBIGIb pm1xBIGb idb id0b id1b id_1b ...
                                        idv1b idm1b iz0b izm_1b izm0b izv0b

                                elseif m~=0 %Non-zonal modified fnALFs                                    
                                    xb=ps1b(:,m);
                                    ixb=ips1b(:,m);
                                    
                                    pmxb=ixb==0;
                                    pm0b=ixb<-1;
                                    pmxBIGIb=(ixb>=-1 & ixb<0);
                                    pmxBIGb=ixb>0;

                                    Pnm(pmxb,1)=xb(pmxb);
                                    Pnm(pm0b,1)=0;
                                    Pnm(pmxBIGIb,1)=xb(pmxBIGIb)*BIGI;
                                    Pnm(pmxBIGb,1)=xb(pmxBIGb)*BIG;
                                   
                                    if m<nmax
                                       yb=xb;
                                       iyb=ixb;

                                       xb=(am(m+1).*t).*yb.*q;
                                       ixb=iyb;
                                       wb=abs(xb);

                                       wBIGSb=wb>=BIGS;
                                       xb(wBIGSb)=xb(wBIGSb)*BIGI;
                                       ixb(wBIGSb)=ixb(wBIGSb)+1;

                                       wBIGSIb=wb<BIGSI;
                                       xb(wBIGSIb)=xb(wBIGSIb)*BIG;
                                       ixb(wBIGSIb)=ixb(wBIGSIb)-1;   

                                       pm1xb=ixb==0;
                                       pm10b=ixb<-1;
                                       pm1xBIGIb=(ixb>=-1 & ixb<0);
                                       pm1xBIGb=ixb>0;

                                       Pnm(pm1xb,2)=xb(pm1xb);
                                       Pnm(pm10b,2)=0;
                                       Pnm(pm1xBIGIb,2)=xb(pm1xBIGIb)*BIGI;
                                       Pnm(pm1xBIGb,2)=xb(pm1xBIGb)*BIG;

                                       for n=m+2:nmax
                                           idb=ixb-iyb;

                                           id0b=idb==0;
                                           id1b=idb==1;
                                           id_1b=idb==-1;
                                           idv1b=idb>1;
                                           idm1b=idb<-1;

                                           zzb(id0b)=(am(n).*t(id0b).*q(id0b)).*xb(id0b)-bm(n).*yb(id0b).*q2(id0b);
                                           izb(id0b,1)=ixb(id0b);

                                           zzb(id1b)=(am(n).*t(id1b).*q(id1b)).*xb(id1b)-bm(n).*(yb(id1b).*q2(id1b).*BIGI);
                                           izb(id1b)=ixb(id1b);

                                           zzb(id_1b)=(am(n).*t(id_1b).*q(id_1b)).*(xb(id_1b).*BIGI)-bm(n).*q2(id_1b).*yb(id_1b);
                                           izb(id_1b)=iyb(id_1b);

                                           zzb(idv1b)=(am(n).*t(idv1b).*q(idv1b)).*xb(idv1b);
                                           izb(idv1b)=ixb(idv1b);

                                           zzb(idm1b)=-bm(n).*yb(idm1b).*q2(idm1b);
                                           izb(idm1b)=iyb(idm1b);

                                           wb=abs(zzb);

                                           wBIGSb=wb>=BIGS;
                                           zzb(wBIGSb)=zzb(wBIGSb)*BIGI;
                                           izb(wBIGSb)=izb(wBIGSb)+1;

                                           wBIGSIb=wb<BIGSI;
                                           zzb(wBIGSIb)=zzb(wBIGSIb)*BIG;
                                           izb(wBIGSIb)=izb(wBIGSIb)-1;

                                           iz0b=izb==0;
                                           izm_1b=izb<-1;
                                           izm0b=(izb>=-1 & izb<0);
                                           izv0b=izb>0;

                                           Pnm(iz0b,n-m+1)=zzb(iz0b);
                                           Pnm(izm_1b,n-m+1)=0;
                                           Pnm(izm0b,n-m+1)=zzb(izm0b)*BIGI;
                                           Pnm(izv0b,n-m+1)=zzb(izv0b)*BIG;                           
                              
                                           yb=xb;
                                           iyb=ixb;
                                           xb=zzb;
                                           ixb=izb;
                                       end   
                                    end
                                end
                            end
                        end

                        %If nmin~=2
                        %======================================================
                        if nmin~=2
                            if LNOFnmin==0
                                if m<2
                                    if por==1
                                        deltaCm(1:(nmin-2))=0;
                                    end

                                    if grav==1
                                        Cm(1:(nmin-2))=0;
                                    end

                                    if geoid==1
                                        HCm(1:(nmin-m))=0;
                                        HSm(1:(nmin-m))=0;
                                    end

                                    if normal==1
                                        if m==0
                                            CElm(1:(nmin-2))=0;
                                        end
                                    end

                                    Sm(1:(nmin-2))=0;
                                else
                                    if por==1
                                        deltaCm(1:(nmin-m))=0;
                                    end

                                    if grav==1
                                        Cm(1:(nmin-m))=0;
                                    end

                                    if geoid==1
                                        HCm(1:(nmin-m))=0;
                                        HSm(1:(nmin-m))=0;
                                    end

                                    Sm(1:(nmin-m))=0;
                                end
                            elseif LNOFnmin==1
                                Pnm(:,1:(nmin-m))=0;
                            end
                        end
                        %======================================================
   
                        %% Computation of the first-order derivatives of modified fnALFs
                        if dALFs==1  
                            if volbaALFs==1 || volbaALFs==3
                                enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                                if m==0 %Zonal modified dALFs
                                    dPnm(:,1)=0.*u;
                                    dPnm(:,2)=sqrt(3)*u.*q;
                                    dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                                elseif m==nmax %Sectorial modified dALFs
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1);
                                else
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                    dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne ALFs
                                end

                            elseif volbaALFs==2
                                enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                                if m==0 %Zonal modified dALFs
                                    dPnm(:,1)=0.*u*1e-280;
                                    dPnm(:,2)=sqrt(3)*u.*q*1e-280;
                                    dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                                elseif m==nmax
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Sectorial modified dALFs
                                else
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                    dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne dALFs
                                end
                            end

                            %Treatment of the dALFs singularity
                            dPnm(singdPnm,:)=0;

                            if ddALFs==1 %If the second-order derivatives of the modified fnALFs are to be computed
                                
                                if m==0 %Zonal modified ddALFs
                                    ddPnm=bsxfun(@times,tu,dPnm)-bsxfun(@times,(0:nmax).*((0:nmax)+1),Pnm);
                                else
                                    ddPnm(:,1:end-m)=bsxfun(@times,tu,dPnm(:,1:end-m))+bsxfun(@times,m^2./u.^2,Pnm(:,1:end-m))-bsxfun(@times,(m:nmax).*((m:nmax)+1),Pnm(:,1:end-m));
                                end
                                                                
                                %Treatment of the ddALFs singularity
                                ddPnm(singddPnm,:)=0;
                            end                                                                               
                        end

                        %% Loop for 1:NF (number of computing functionals)
                        for i=1:pocetpar                   
                            if volbapar(i)==1       
                            elseif volbapar(i)==2 %Deflection of the vertical eta                       
                                
                                if m<2
                                    Lm=m*Pnm(:,(3-m):(end-m));
                                else
                                    Lm=m*Pnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3 %Computation using the standard 
                                    %forward column method or extended-range arithmetic
                                    if m==nmax
                                        Aeta=zeros(length(fiG),nmax+1);
                                        Beta=zeros(length(fiG),nmax+1);
                                    end

                                    Aeta(:,m+1)=Lm*deltaCm;
                                    Beta(:,m+1)=Lm*Sm;                        
                                elseif volbaALFs==2 %Computation using the modified forward column method combined with Horner's scheme
                                    if m==nmax
                                        eta=zeros(length(fiG),length(lambda));
                                    end

                                    eta=bsxfun(@times,eta,u)+(-Lm*deltaCm*sin(m*lambda')+Lm*Sm*cos(m*lambda'));
                                end

                            elseif volbapar(i)==3 %Deflection of the vertical xi

                                if m==nmax
                                    Aksi=zeros(length(fiG),nmax+1);
                                    Bksi=zeros(length(fiG),nmax+1);
                                end

                                if m<2
                                    dLm=dPnm(:,(3-m):(end-m));
                                else
                                    dLm=dPnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        Aksi=zeros(length(fiG),nmax+1);
                                        Bksi=zeros(length(fiG),nmax+1);
                                    end

                                    Aksi(:,m+1)=dLm*deltaCm;
                                    Bksi(:,m+1)=dLm*Sm; 
                                elseif volbaALFs==2
                                    if m==nmax
                                        ksi=zeros(length(fiG),length(lambda));
                                    end

                                    ksi=bsxfun(@times,ksi,u)+(dLm*deltaCm*cos(m*lambda')+dLm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==4 %Deflection of the vertical Theta

                                if m<2
                                    Lm=m*Pnm(:,(3-m):(end-m));
                                    dLm=dPnm(:,(3-m):(end-m));
                                else
                                    Lm=m*Pnm(:,1:(nmax-m+1));
                                    dLm=dPnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        ATeta=zeros(length(fiG),nmax+1);
                                        BTeta=zeros(length(fiG),nmax+1);
                                        ATksi=ATeta;
                                        BTksi=BTeta;
                                    end

                                    ATeta(:,m+1)=Lm*deltaCm;
                                    BTeta(:,m+1)=Lm*Sm;
                                    ATksi(:,m+1)=dLm*deltaCm;
                                    BTksi(:,m+1)=dLm*Sm; 
                                elseif volbaALFs==2
                                    if m==nmax
                                        Teta=zeros(length(fiG),length(lambda));
                                        Tksi=zeros(length(fiG),length(lambda));
                                    end

                                    Teta=bsxfun(@times,Teta,u)+(-Lm*deltaCm*sin(m*lambda')+Lm*Sm*cos(m*lambda'));
                                    Tksi=bsxfun(@times,Tksi,u)+(dLm*deltaCm*cos(m*lambda')+dLm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==5 %Disturbing potential

                                if m<2
                                    Lm=Pnm(:,(3-m):(end-m));
                                else
                                    Lm=Pnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AT=zeros(length(fiG),nmax+1);
                                        BT=zeros(length(fiG),nmax+1);
                                    end

                                    AT(:,m+1)=Lm*deltaCm;
                                    BT(:,m+1)=Lm*Sm;                                    
                                elseif volbaALFs==2
                                    if m==nmax
                                        T=zeros(length(fiG),length(lambda));
                                    end

                                    T=bsxfun(@times,T,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                                
                                if m==nmax                               
                                    ampl_Trr=((2:nmax)+1).*((2:nmax)+2);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_Trr);
                                    Lmff=ddPnm(:,(3-m):(end-m));
                                    Lmll=m^2*Pnm(:,(3-m):(end-m));
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Trr(:,(m-1):end));
                                    Lmff=ddPnm(:,1:(nmax-m+1));
                                    Lmll=m^2*Pnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        ATrr=zeros(length(fiG),nmax+1);
                                        BTrr=zeros(length(fiG),nmax+1);
                                        ATff=ATrr;
                                        BTff=ATrr;
                                        ATll=ATrr;
                                        BTll=ATrr;
                                    end

                                    ATrr(:,m+1)=Lm*deltaCm;
                                    BTrr(:,m+1)=Lm*Sm;  
                                    
                                    ATff(:,m+1)=Lmff*deltaCm;
                                    BTff(:,m+1)=Lmff*Sm;
                                    
                                    ATll(:,m+1)=Lmll*deltaCm;
                                    BTll(:,m+1)=Lmll*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Trr=zeros(length(fiG),length(lambda));
                                        Tff=Trr;
                                        Tll=Trr;
                                    end

                                    Trr=bsxfun(@times,Trr,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                    Tff=bsxfun(@times,Tff,u)+(Lmff*deltaCm*cos(m*lambda')+Lmff*Sm*sin(m*lambda'));
                                    Tll=bsxfun(@times,Tll,u)+(Lmll*deltaCm*cos(m*lambda')+Lmll*Sm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                                
                                if m==nmax                               
                                    ampl_Tfl=(2:nmax)+1;
                                end
                                
                                if m<2
                                    Lmrf=bsxfun(@times,dPnm(:,(3-m):(end-m)),ampl_Tfl);
                                    Lmrl=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),ampl_Tfl);
                                    Lmfl=m*dPnm(:,(3-m):(end-m));
                                else
                                    Lmrf=bsxfun(@times,dPnm(:,1:(nmax-m+1)),ampl_Tfl(:,(m-1):end));
                                    Lmrl=bsxfun(@times,m*Pnm(:,1:(nmax-m+1)),ampl_Tfl(:,(m-1):end));
                                    Lmfl=m*dPnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        ATrf=zeros(length(fiG),nmax+1);
                                        BTrf=zeros(length(fiG),nmax+1);
                                        ATrl=ATrf;
                                        BTrl=ATrf;
                                        ATfl=ATrf;
                                        BTfl=ATrf;
                                    end

                                    ATrf(:,m+1)=Lmrf*deltaCm;
                                    BTrf(:,m+1)=Lmrf*Sm;  
                                    
                                    ATrl(:,m+1)=Lmrl*deltaCm;
                                    BTrl(:,m+1)=Lmrl*Sm;
                                    
                                    ATfl(:,m+1)=Lmfl*deltaCm;
                                    BTfl(:,m+1)=Lmfl*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Trf=zeros(length(fiG),length(lambda));
                                        Trl=Trf;
                                        Tfl=Trf;
                                    end

                                    Trf=bsxfun(@times,Trf,u)+(Lmrf*deltaCm*cos(m*lambda')+Lmrf*Sm*sin(m*lambda'));
                                    Trl=bsxfun(@times,Trl,u)+(Lmrl*Sm*cos(m*lambda')-Lmrl*deltaCm*sin(m*lambda'));
                                    Tfl=bsxfun(@times,Tfl,u)+(Lmfl*Sm*cos(m*lambda')-Lmfl*deltaCm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                                
                                if m==nmax                               
                                    ampl_Tzz=((2:nmax)+1).*((2:nmax)+2);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_Tzz);
                                       
                                    %Txx
                                    bnm=((2:nmax)+m+1).*((2:nmax)+m+2)./2./(m+1);
                                    LmTxx1=bsxfun(@times,Pnm(:,(3-m):(end-m)),bnm-((2:nmax)+1).*((2:nmax)+2));
                                    LmTyy1=bsxfun(@times,Pnm(:,(3-m):(end-m)),bnm);
                                    
                                    if m==0
                                        anm=sqrt(2)/4.*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)-(m+2)+2);
                                    else
                                        anm=1./4.*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)-(m+2)+2);
                                    end
                                    
                                    LmTxx2=bsxfun(@times,Pnm(:,(3-m):(end-m)),anm);
                                    LmTxx3=zeros(length(fi),nmax-1); %Coefficients cnm are equal to zeros, if m==0 a m==1
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Tzz(:,(m-1):end));
                                 
                                    %Txx
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    LmTxx1=bsxfun(@times,Pnm(:,1:(nmax-m+1)),bnm-((m:nmax)+1).*((m:nmax)+2));
                                    LmTyy1=bsxfun(@times,Pnm(:,1:(nmax-m+1)),bnm);

                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    LmTxx2=bsxfun(@times,Pnm(:,1:(nmax-m+1)),anm);
                                       
                                    if m==2
                                        cnm=sqrt(2)/4.*sqrt((2:nmax).^2-(m-2+1).^2).*sqrt((2:nmax)-(m-2)).*sqrt((2:nmax)+m-2+2);
                                    else
                                        cnm=1./4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                    end
                                    
                                    LmTxx3=bsxfun(@times,Pnm(:,1:(nmax-m+1)),cnm);                                    
                                end
                   
                                if m==nmax
                                    ATzz=zeros(length(fiG),nmax+1);
                                    BTzz=zeros(length(fiG),nmax+1);
                                    ATxx1=ATzz;
                                    BTxx1=ATzz;
                                    ATxx2=ATzz;
                                    BTxx2=ATzz;
                                    ATxx3=ATzz;
                                    BTxx3=BTzz;
                                    ATyy1=ATzz;
                                    BTyy1=BTzz;
                                end

                                ATzz(:,m+1)=Lm*deltaCm;
                                BTzz(:,m+1)=Lm*Sm;  
                                    
                                %Txx
                                ATxx1(:,m+1)=LmTxx1*deltaCm;
                                BTxx1(:,m+1)=LmTxx1*Sm;
 
                                if radenie==0
                                    if m<2
                                        ATxx2(:,m+1)=LmTxx2*deltaC(index+m+2);
                                        BTxx2(:,m+1)=LmTxx2*S(index+m+2);

                                        ATxx3(:,m+1)=0;
                                        BTxx3(:,m+1)=0;
                                    elseif m>nmax-2
                                        ATxx2(:,m+1)=0;
                                        BTxx2(:,m+1)=0;
                                           
                                        ATxx3(:,m+1)=LmTxx3*deltaC(index((m-1):end)+m-2);
                                        BTxx3(:,m+1)=LmTxx3*S(index((m-1):end)+m-2);
                                    else
                                        ATxx2(:,m+1)=LmTxx2*deltaC(index((m-1):end)+m+2);
                                        BTxx2(:,m+1)=LmTxx2*S(index((m-1):end)+m+2);

                                        ATxx3(:,m+1)=LmTxx3*deltaC(index((m-1):end)+m-2);
                                        BTxx3(:,m+1)=LmTxx3*S(index((m-1):end)+m-2);
                                    end
                                elseif radenie==1
                                    if m<2
                                        ATxx2(:,m+1)=LmTxx2*[zeros(m,1);deltaC((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))];
                                        BTxx2(:,m+1)=LmTxx2*[zeros(m,1);S((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))];

                                        ATxx3(:,m+1)=0;
                                        BTxx3(:,m+1)=0;
                                    elseif m>nmax-2
                                        ATxx2(:,m+1)=0;
                                        BTxx2(:,m+1)=0;
       
                                        if m==2
                                            ATxx3(:,m+1)=LmTxx3*deltaC(1:(nmax-1));
                                            BTxx3(:,m+1)=0;
                                        else
                                            ATxx3(:,m+1)=LmTxx3*deltaC((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            BTxx3(:,m+1)=LmTxx3*S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3))); 
                                        end                                                                                 
                                    else                                        
                                        ATxx2(:,m+1)=LmTxx2*[0;0;deltaC((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                        BTxx2(:,m+1)=LmTxx2*[0;0;S((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];

                                        if m==2
                                            ATxx3(:,m+1)=LmTxx3*deltaC(1:(nmax-1));
                                            BTxx3(:,m+1)=0;
                                        else
                                            ATxx3(:,m+1)=LmTxx3*deltaC((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            BTxx3(:,m+1)=LmTxx3*S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                    end
                                end

                                %Tyy
                                ATyy1(:,m+1)=LmTyy1*deltaCm;
                                BTyy1(:,m+1)=LmTyy1*Sm;
                                    
                                %Modified forward column method cannot be
                                %applied 
                                
                            elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                                                                
                                if m<2                                       
                                    if m==0
                                        betanm=((2:nmax)+2)./2.*sqrt(1+ones(1,nmax-1)).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)-(m+1)+1);
                                        gamanm=zeros(1,nmax-1);
                                        
                                        %Txy
                                        dnm=1/4.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt(2).*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)+m+2-2);                                      
                                        LmTxy1=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),dnm);

                                        LmTxy2=zeros(length(fi),nmax-1);
                                        
                                        LmTxy3=zeros(length(fi),nmax-1);
                                        
                                        %Tyz
                                        minm=((2:nmax)+2)./2.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt(2).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)+m+1-1);
                                        LmTyz1=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),minm);
                                        
                                        LmTyz2=zeros(length(fi),nmax-1);
                                    else
                                        betanm=((2:nmax)+2)./2.*sqrt((2:nmax)+m+1).*sqrt((2:nmax)-(m+1)+1);
                                        gamanm=-((2:nmax)+2).*sqrt((2:nmax).*((2:nmax)+1)./2);
                                        
                                        %Txy
                                        dnm=1/4.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)+m+2-2);                                        
                                        LmTxy1=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),dnm);
                                        
                                        gnm=-1/4*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax)+1).*sqrt((2:nmax)-1).*((2:nmax)+2);
                                        LmTxy2=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),gnm);
                                        
                                        LmTxy3=zeros(length(fi),nmax-1);
                                        
                                        %Tyz
                                        minm=((2:nmax)+2)./2.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)+m+1-1);
                                        LmTyz1=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),minm);
                                        
                                        LmTyz2=zeros(length(fi),nmax-1);
                                    end
                                    
                                    %Txz
                                    LmTxz1=bsxfun(@times,Pnm(:,(3-m):(end-m)),betanm);
                                    LmTxz2=bsxfun(@times,Pnm(:,(3-m):(end-m)),gamanm); 
                                else
                                    %Txz
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                                                                                                                               
                                    LmTxz1=bsxfun(@times,Pnm(:,1:(nmax-m+1)),betanm);
                                    LmTxz2=bsxfun(@times,Pnm(:,1:(nmax-m+1)),gamanm);
                                                                     
                                    %Txy
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    LmTxy1=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],dnm);
                                    
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    LmTxy2=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],gnm);
                                    
                                    if m==2
                                        %Txy
                                        LmTxy3=zeros(length(fi),nmax-1);    
                                    elseif m==3
                                        %Txy
                                        hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-3).*sqrt((m:nmax)-2).*sqrt((m:nmax)-1).*sqrt((m:nmax)+2);
                                        LmTxy3=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],hnm);
                                    else
                                        %Txy
                                        hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)-(m-2)-2);
                                        LmTxy3=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],hnm);
                                    end
                                    
                                    %Tyz
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    LmTyz1=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],minm);

                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                    LmTyz2=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],ninm);
                                end

                                if m==nmax
                                    ATxz1=zeros(length(fiG),nmax+1);
                                    BTxz1=zeros(length(fiG),nmax+1);
                                    ATxz2=ATxz1;
                                    BTxz2=ATxz1;
                                    ATxy1=ATxz1;
                                    BTxy1=ATxz1;
                                    ATxy2=ATxz1;
                                    BTxy2=ATxz1;
                                    ATxy3=ATxz1;
                                    BTxy3=ATxz1;
                                    ATyz1=ATxz1;
                                    BTyz1=ATxz1;
                                    ATyz2=ATxz1;
                                    BTyz2=ATxz1;
                                end

                                ATxy2(:,m+1)=(LmTxy2*deltaCm).*q;
                                BTxy2(:,m+1)=(LmTxy2*Sm).*q;
                                    
                                if radenie==0
                                    if m<2
                                        ATxz1(:,m+1)=LmTxz1*deltaC(index+m+1);
                                        BTxz1(:,m+1)=LmTxz1*S(index+m+1);

                                        if m==1
                                            %Txz
                                            ATxz2(:,m+1)=LmTxz2*deltaC(index+m-1);
                                            BTxz2(:,m+1)=LmTxz2*S(index+m-1);
                                                
                                            %Tyz
                                            ATyz2(:,m+1)=LmTyz2*deltaC(index+m-1).*q;
                                            BTyz2(:,m+1)=LmTyz2*S(index+m-1).*q;
                                        else
                                            %Txz
                                            ATxz2(:,m+1)=0;
                                            BTxz2(:,m+1)=0;
                                                
                                            %Tyz
                                            ATyz2(:,m+1)=0;
                                            BTyz2(:,m+1)=0;
                                        end
                                        
                                        %Txy
                                        ATxy1(:,m+1)=(LmTxy1*deltaC(index+m+2)).*q;  
                                        BTxy1(:,m+1)=(LmTxy1*S(index+m+2)).*q; 
                                           
                                        ATxy3(:,m+1)=LmTxy3*zeros(nmax-2+1,1);
                                        BTxy3(:,m+1)=LmTxy3*zeros(nmax-2+1,1);
                                            
                                        %Tyz
                                        ATyz1(:,m+1)=(LmTyz1*deltaC(index+m+1)).*q;  
                                        BTyz1(:,m+1)=(LmTyz1*S(index+m+1)).*q; 
                                    elseif m>nmax-2
                                            
                                        if m==nmax
                                            %Txz
                                            ATxz1(:,m+1)=0;
                                            BTxz1(:,m+1)=0;
                                                
                                            %Tyz
                                            ATyz1(:,m+1)=0;
                                            BTyz1(:,m+1)=0;
                                        else
                                            %Txz
                                            ATxz1(:,m+1)=LmTxz1*deltaC(index((m-1):end)+m+1);
                                            BTxz1(:,m+1)=LmTxz1*S(index((m-1):end)+m+1);
                                                
                                            %Tyz
                                            ATyz1(:,m+1)=LmTyz1*deltaC(index((m-1):end)+m+1).*q;
                                            BTyz1(:,m+1)=LmTyz1*S(index((m-1):end)+m+1).*q;
                                        end
                                           
                                        ATxz2(:,m+1)=LmTxz2*deltaC(index((m-1):end)+m-1);
                                        BTxz2(:,m+1)=LmTxz2*S(index((m-1):end)+m-1);

                                        %Txy
                                        ATxy1(:,m+1)=0;
                                        BTxy1(:,m+1)=0;
                                            
                                        ATxy3(:,m+1)=(LmTxy3*deltaC(index((m-1):end)+m-2)).*q;
                                        BTxy3(:,m+1)=(LmTxy3*S(index((m-1):end)+m-2)).*q;  
                                            
                                        %Tyz
                                        ATyz2(:,m+1)=LmTyz2*deltaC(index((m-1):end)+m-1).*q;
                                        BTyz2(:,m+1)=LmTyz2*S(index((m-1):end)+m-1).*q;
                                    else
                                        ATxz1(:,m+1)=LmTxz1*deltaC(index((m-1):end)+m+1);
                                        BTxz1(:,m+1)=LmTxz1*S(index((m-1):end)+m+1);

                                        ATxz2(:,m+1)=LmTxz2*deltaC(index((m-1):end)+m-1);
                                        BTxz2(:,m+1)=LmTxz2*S(index((m-1):end)+m-1);

                                        %Txy
                                        ATxy1(:,m+1)=(LmTxy1*deltaC(index((m-1):end)+m+2)).*q;
                                        BTxy1(:,m+1)=(LmTxy1*S(index((m-1):end)+m+2)).*q;
                                            
                                        ATxy3(:,m+1)=(LmTxy3*deltaC(index((m-1):end)+m-2)).*q;
                                        BTxy3(:,m+1)=(LmTxy3*S(index((m-1):end)+m-2)).*q;
                                            
                                        %Tyz
                                        ATyz1(:,m+1)=LmTyz1*deltaC(index((m-1):end)+m+1).*q;
                                        BTyz1(:,m+1)=LmTyz1*S(index((m-1):end)+m+1).*q;
                                           
                                        ATyz2(:,m+1)=LmTyz2*deltaC(index((m-1):end)+m-1).*q;
                                        BTyz2(:,m+1)=LmTyz2*S(index((m-1):end)+m-1).*q;
                                    end
                                elseif radenie==1
                                    if m<2
                                        %Txy
                                        ATxy1(:,m+1)=(LmTxy1*[zeros(m,1);deltaC((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))]).*q;
                                        BTxy1(:,m+1)=(LmTxy1*[zeros(m,1);S((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))]).*q;

                                        ATxy3(:,m+1)=0;
                                        BTxy3(:,m+1)=0;
                                            
                                        %Txz
                                        ATxz1(:,m+1)=(LmTxz1*deltaC((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)))); 
                                        BTxz1(:,m+1)=(LmTxz1*S((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)))); 
                                            
                                        %Tyz
                                        ATyz1(:,m+1)=(LmTyz1*deltaC((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)))).*q; 
                                        BTyz1(:,m+1)=(LmTyz1*S((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)))).*q; 
                                            
                                        if m==1
                                            %Txz
                                            ATxz2(:,m+1)=LmTxz2*deltaC(1:(nmax-1));
                                            BTxz2(:,m+1)=0;
                                                
                                            %Tyz
                                            ATyz2(:,m+1)=(LmTyz2*deltaC(1:(nmax-1))).*q;
                                            BTyz2(:,m+1)=0;
                                        else
                                            %Txz
                                            ATxz2(:,m+1)=0;
                                            BTxz2(:,m+1)=0;
                                                
                                            %Tyz
                                            ATyz2(:,m+1)=0;
                                            BTyz2(:,m+1)=0;
                                        end
                                    elseif m>nmax-2
                                        if m==2
                                            %Txy
                                            ATxy3(:,m+1)=(LmTxy3*deltaC(1:(nmax-1))).*q;
                                            BTxy3(:,m+1)=0;
                                        else
                                            %Txy
                                            ATxy3(:,m+1)=(LmTxy3*deltaC((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)))).*q;
                                            BTxy3(:,m+1)=(LmTxy3*S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)))).*q; 
                                        end   
                                            
                                        ATxy1(:,m+1)=0;
                                        BTxy1(:,m+1)=0;
                                            
                                        if m==nmax
                                            %Txz
                                            ATxz1(:,m+1)=0;    
                                            BTxz1(:,m+1)=0;
                                                
                                            %Tyz
                                            ATyz1(:,m+1)=0;    
                                            BTyz1(:,m+1)=0;
                                        else
                                            %Txz
                                            ATxz1(:,m+1)=LmTxz1*[0;deltaC((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                            BTxz1(:,m+1)=LmTxz1*[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                                
                                            %Tyz
                                            ATyz1(:,m+1)=(LmTyz1*[0;deltaC((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]).*q;
                                            BTyz1(:,m+1)=(LmTyz1*[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]).*q;
                                        end

                                        ATxz2(:,m+1)=LmTxz2*deltaC((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                        BTxz2(:,m+1)=LmTxz2*S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                          
                                        %Tyz
                                        ATyz2(:,m+1)=LmTyz2*deltaC((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))).*q;
                                        BTyz2(:,m+1)=LmTyz2*S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))).*q;
                                    else                   
                                        %Txy
                                        ATxy1(:,m+1)=(LmTxy1*[0;0;deltaC((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))]).*q;
                                        BTxy1(:,m+1)=(LmTxy1*[0;0;S((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))]).*q;

                                        if m==2
                                            ATxy3(:,m+1)=(LmTxy3*deltaC(1:(nmax-1))).*q;
                                            BTxy3(:,m+1)=0;
                                        else
                                            ATxy3(:,m+1)=(LmTxy3*deltaC((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)))).*q;
                                            BTxy3(:,m+1)=(LmTxy3*S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)))).*q;
                                        end
                                        
                                        %Txz
                                        ATxz1(:,m+1)=(LmTxz1*[0;deltaC((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]);
                                        BTxz1(:,m+1)=(LmTxz1*[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]);
                                            
                                        ATxz2(:,m+1)=LmTxz2*deltaC((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                        BTxz2(:,m+1)=LmTxz2*S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                            
                                        %Tyz
                                        ATyz2(:,m+1)=LmTyz2*deltaC((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))).*q;
                                        BTyz2(:,m+1)=LmTyz2*S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))).*q;
                                            
                                        ATyz1(:,m+1)=(LmTyz1*[0;deltaC((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]).*q;
                                        BTyz1(:,m+1)=(LmTyz1*[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]).*q;
                                    end
                                end
                                
                                %Modified forward column method cannot be
                                %applied
                                
                            elseif volbapar(i)==10 %Geoid undulation

                                if m==nmax
                                    amplH=zeros(length(fiG),nmax+1); %Damping factor
                                    for n=0:nmax
                                        amplH(:,n+1)=1./((R./r).^n);

                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFFIED
                                        % fnALFS
                                    end
                                end
                                        
                                if m<2
                                    Lm=Pnm(:,(3-m):(end-m));
                                    LmH=Pnm(:,1:(nmax-m+1)).*amplH(:,(m+1):end);
                                else
                                    Lm=Pnm(:,1:(nmax-m+1));
                                    LmH=Pnm(:,1:(nmax-m+1)).*amplH(:,(m+1):end);
                                end
                                
                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AN1c=zeros(length(fiG),nmax);
                                        BN1c=AN1c;
                                        AH=AN1c;
                                        BH=AN1c;
                                    end  
                                    
                                    AN1c(:,m+1)=Lm*deltaCm;
                                    BN1c(:,m+1)=Lm*Sm;  
                                    AH(:,m+1)=LmH*HCm;
                                    BH(:,m+1)=LmH*HSm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        N1c=zeros(length(fiG),length(lambda));
                                        H=N1c;
                                    end

                                    N1c=bsxfun(@times,N1c,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                    H=bsxfun(@times,H,u)+(LmH*HCm*cos(m*lambda')+LmH*HSm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==11 %Gravitational potential

                                if m<2
                                    Lm=Pnm(:,(3-m):(end-m));
                                else
                                    Lm=Pnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AV=zeros(length(fiG),nmax+1);
                                        BV=zeros(length(fiG),nmax+1);
                                    end

                                    AV(:,m+1)=Lm*Cm;
                                    BV(:,m+1)=Lm*Sm;  
                                elseif volbaALFs==2
                                    if m==nmax
                                        V=zeros(length(fiG),length(lambda));
                                    end

                                    V=bsxfun(@times,V,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                                
                                if m==nmax                               
                                    ampl_Vrr=((2:nmax)+1).*((2:nmax)+2);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_Vrr);
                                    Lmff=ddPnm(:,(3-m):(end-m));
                                    Lmll=m^2*Pnm(:,(3-m):(end-m));
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vrr(:,(m-1):end));
                                    Lmff=ddPnm(:,1:(nmax-m+1));
                                    Lmll=m^2*Pnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AVrr=zeros(length(fiG),nmax+1);
                                        BVrr=zeros(length(fiG),nmax+1);
                                        AVff=AVrr;
                                        BVff=AVrr;
                                        AVll=AVrr;
                                        BVll=AVrr;
                                    end

                                    AVrr(:,m+1)=Lm*Cm;
                                    BVrr(:,m+1)=Lm*Sm;  
                                    
                                    AVff(:,m+1)=Lmff*Cm;
                                    BVff(:,m+1)=Lmff*Sm;
                                    
                                    AVll(:,m+1)=Lmll*Cm;
                                    BVll(:,m+1)=Lmll*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Vrr=zeros(length(fiG),length(lambda));
                                        Vff=Vrr;
                                        Vll=Vrr;
                                    end

                                    Vrr=bsxfun(@times,Vrr,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                    Vff=bsxfun(@times,Vff,u)+(Lmff*Cm*cos(m*lambda')+Lmff*Sm*sin(m*lambda'));
                                    Vll=bsxfun(@times,Vll,u)+(Lmll*Cm*cos(m*lambda')+Lmll*Sm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                                
                                if m==nmax                               
                                    ampl_Vfl=((2:nmax)+1);
                                end
                                
                                if m<2
                                    Lmrf=bsxfun(@times,dPnm(:,(3-m):(end-m)),ampl_Vfl);
                                    Lmrl=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),ampl_Vfl);
                                    Lmfl=m*dPnm(:,(3-m):(end-m));
                                else
                                    Lmrf=bsxfun(@times,dPnm(:,1:(nmax-m+1)),ampl_Vfl(:,(m-1):end));
                                    Lmrl=bsxfun(@times,m*Pnm(:,1:(nmax-m+1)),ampl_Vfl(:,(m-1):end));
                                    Lmfl=m*dPnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AVrf=zeros(length(fiG),nmax+1);
                                        BVrf=zeros(length(fiG),nmax+1);
                                        AVrl=AVrf;
                                        BVrl=AVrf;
                                        AVfl=AVrf;
                                        BVfl=AVrf;
                                    end

                                    AVrf(:,m+1)=Lmrf*Cm;
                                    BVrf(:,m+1)=Lmrf*Sm;  
                                    
                                    AVrl(:,m+1)=Lmrl*Cm;
                                    BVrl(:,m+1)=Lmrl*Sm;
                                    
                                    AVfl(:,m+1)=Lmfl*Cm;
                                    BVfl(:,m+1)=Lmfl*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Vrf=zeros(length(fiG),length(lambda));
                                        Vrl=Vrf;
                                        Vfl=Vrf;
                                    end

                                    Vrf=bsxfun(@times,Vrf,u)+(Lmrf*Cm*cos(m*lambda')+Lmrf*Sm*sin(m*lambda'));
                                    Vrl=bsxfun(@times,Vrl,u)+(Lmrl*Sm*cos(m*lambda')-Lmrl*Cm*sin(m*lambda'));
                                    Vfl=bsxfun(@times,Vfl,u)+(Lmfl*Sm*cos(m*lambda')-Lmfl*Cm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz

                                if m==nmax                               
                                    ampl_Vzz=((2:nmax)+1).*((2:nmax)+2);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_Vzz);
                                       
                                    %Vxx
                                    bnm=((2:nmax)+m+1).*((2:nmax)+m+2)./2./(m+1);
                                    LmVxx1=bsxfun(@times,Pnm(:,(3-m):(end-m)),bnm-((2:nmax)+1).*((2:nmax)+2));
                                    LmVyy1=bsxfun(@times,Pnm(:,(3-m):(end-m)),bnm);
                                    
                                    if m==0
                                        anm=sqrt(2)/4.*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)-(m+2)+2);
                                    else
                                        anm=1./4.*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)-(m+2)+2);
                                    end
                                    
                                    LmVxx2=bsxfun(@times,Pnm(:,(3-m):(end-m)),anm);
                                    LmVxx3=zeros(length(fi),nmax-1); %Coefficients cnm are equal to zero, if m==0 a m==1
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Vzz(:,(m-1):end));
                                 
                                    %Vxx
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    LmVxx1=bsxfun(@times,Pnm(:,1:(nmax-m+1)),bnm-((m:nmax)+1).*((m:nmax)+2));
                                    LmVyy1=bsxfun(@times,Pnm(:,1:(nmax-m+1)),bnm);

                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    LmVxx2=bsxfun(@times,Pnm(:,1:(nmax-m+1)),anm);
                                       
                                    if m==2
                                        cnm=sqrt(2)/4.*sqrt((2:nmax).^2-(m-2+1).^2).*sqrt((2:nmax)-(m-2)).*sqrt((2:nmax)+m-2+2);
                                    else
                                        cnm=1./4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                    end
                                    
                                    LmVxx3=bsxfun(@times,Pnm(:,1:(nmax-m+1)),cnm);                                    
                                end
                   
                                if m==nmax
                                    AVzz=zeros(length(fiG),nmax+1);
                                    BVzz=zeros(length(fiG),nmax+1);
                                    AVxx1=AVzz;
                                    BVxx1=AVzz;
                                    AVxx2=AVzz;
                                    BVxx2=AVzz;
                                    AVxx3=AVzz;
                                    BVxx3=BVzz;
                                    AVyy1=AVzz;
                                    BVyy1=BVzz;
                                end

                                AVzz(:,m+1)=Lm*Cm;
                                BVzz(:,m+1)=Lm*Sm;  
                                    
                                %Vxx
                                AVxx1(:,m+1)=LmVxx1*Cm;
                                BVxx1(:,m+1)=LmVxx1*Sm;
 
                                if radenie==0
                                    if m<2
                                        AVxx2(:,m+1)=LmVxx2*C(index+m+2);
                                        BVxx2(:,m+1)=LmVxx2*S(index+m+2);

                                        AVxx3(:,m+1)=0;
                                        BVxx3(:,m+1)=0;
                                    elseif m>nmax-2
                                        AVxx2(:,m+1)=0;
                                        BVxx2(:,m+1)=0;
                                           
                                        AVxx3(:,m+1)=LmVxx3*C(index((m-1):end)+m-2);
                                        BVxx3(:,m+1)=LmVxx3*S(index((m-1):end)+m-2);
                                    else
                                        AVxx2(:,m+1)=LmVxx2*C(index((m-1):end)+m+2);
                                        BVxx2(:,m+1)=LmVxx2*S(index((m-1):end)+m+2);

                                        AVxx3(:,m+1)=LmVxx3*C(index((m-1):end)+m-2);
                                        BVxx3(:,m+1)=LmVxx3*S(index((m-1):end)+m-2);
                                    end
                                elseif radenie==1
                                    if m<2
                                        AVxx2(:,m+1)=LmVxx2*[zeros(m,1);C((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))];
                                        BVxx2(:,m+1)=LmVxx2*[zeros(m,1);S((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))];

                                        AVxx3(:,m+1)=0;
                                        BVxx3(:,m+1)=0;
                                    elseif m>nmax-2
                                        AVxx2(:,m+1)=0;
                                        BVxx2(:,m+1)=0;
       
                                        if m==2
                                            AVxx3(:,m+1)=LmVxx3*C(1:(nmax-1));
                                            BVxx3(:,m+1)=0;
                                        else
                                            AVxx3(:,m+1)=LmVxx3*C((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            BVxx3(:,m+1)=LmVxx3*S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3))); 
                                        end                                                                                 
                                    else                                        
                                        AVxx2(:,m+1)=LmVxx2*[0;0;C((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                        BVxx2(:,m+1)=LmVxx2*[0;0;S((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];

                                        if m==2
                                            AVxx3(:,m+1)=LmVxx3*C(1:(nmax-1));
                                            BVxx3(:,m+1)=0;
                                        else
                                            AVxx3(:,m+1)=LmVxx3*C((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            BVxx3(:,m+1)=LmVxx3*S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                    end
                                end

                                %Vyy
                                AVyy1(:,m+1)=LmVyy1*Cm;
                                BVyy1(:,m+1)=LmVyy1*Sm;
                                    
                                %Modified forward column method cannot be
                                %applied
                                
                            elseif volbapar(i)==15 %Gravitational tensor Vxy_Vxz_Vyz
                                
                                if m<2                                       
                                    if m==0
                                        betanm=((2:nmax)+2)./2.*sqrt(1+ones(1,nmax-1)).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)-(m+1)+1);
                                        gamanm=zeros(1,nmax-1);
                                        
                                        %Vxy
                                        dnm=1/4.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt(2).*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)+m+2-2);                                      
                                        LmVxy1=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),dnm);

                                        LmVxy2=zeros(length(fi),nmax-1);
                                        
                                        LmVxy3=zeros(length(fi),nmax-1);
                                        
                                        %Vyz
                                        minm=((2:nmax)+2)./2.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt(2).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)+m+1-1);
                                        LmVyz1=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),minm);
                                        
                                        LmVyz2=zeros(length(fi),nmax-1);
                                    else
                                        betanm=((2:nmax)+2)./2.*sqrt((2:nmax)+m+1).*sqrt((2:nmax)-(m+1)+1);
                                        gamanm=-((2:nmax)+2).*sqrt((2:nmax).*((2:nmax)+1)./2);
                                        
                                        %Vxy
                                        dnm=1/4.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)+m+2-2);                                        
                                        LmVxy1=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),dnm);
                                        
                                        gnm=-1/4*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax)+1).*sqrt((2:nmax)-1).*((2:nmax)+2);
                                        LmVxy2=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),gnm);
                                        
                                        LmVxy3=zeros(length(fi),nmax-1);
                                        
                                        %Vyz
                                        minm=((2:nmax)+2)./2.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)+m+1-1);
                                        LmVyz1=bsxfun(@times,Pnm(:,(3-m-1):(end-m-1)),minm);
                                        
                                        LmVyz2=zeros(length(fi),nmax-1);
                                    end
                                    
                                    %Vxz
                                    LmVxz1=bsxfun(@times,Pnm(:,(3-m):(end-m)),betanm);
                                    LmVxz2=bsxfun(@times,Pnm(:,(3-m):(end-m)),gamanm); 
                                else
                                    %Vxz
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                                                                                                                               
                                    LmVxz1=bsxfun(@times,Pnm(:,1:(nmax-m+1)),betanm);
                                    LmVxz2=bsxfun(@times,Pnm(:,1:(nmax-m+1)),gamanm);
                                                                     
                                    %Vxy
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    LmVxy1=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],dnm);
                                    
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    LmVxy2=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],gnm);
                                    
                                    if m==2
                                        %Vxy
                                        LmVxy3=zeros(length(fi),nmax-1);    
                                    elseif m==3
                                        %Vxy
                                        hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-3).*sqrt((m:nmax)-2).*sqrt((m:nmax)-1).*sqrt((m:nmax)+2);
                                        LmVxy3=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],hnm);
                                    else
                                        %Vxy
                                        hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)-(m-2)-2);
                                        LmVxy3=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],hnm);
                                    end
                                    
                                    %Vyz
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    LmVyz1=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],minm);

                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                    LmVyz2=bsxfun(@times,[zeros(length(fiG),1) Pnm(:,1:(nmax-m))],ninm);
                                end

                                if m==nmax
                                    AVxz1=zeros(length(fiG),nmax+1);
                                    BVxz1=zeros(length(fiG),nmax+1);
                                    AVxz2=AVxz1;
                                    BVxz2=AVxz1;
                                    AVxy1=AVxz1;
                                    BVxy1=AVxz1;
                                    AVxy2=AVxz1;
                                    BVxy2=AVxz1;
                                    AVxy3=AVxz1;
                                    BVxy3=AVxz1;
                                    AVyz1=AVxz1;
                                    BVyz1=AVxz1;
                                    AVyz2=AVxz1;
                                    BVyz2=AVxz1;
                                end

                                AVxy2(:,m+1)=(LmVxy2*Cm).*q;
                                BVxy2(:,m+1)=(LmVxy2*Sm).*q;
                                    
                                if radenie==0
                                    if m<2
                                        AVxz1(:,m+1)=LmVxz1*C(index+m+1);
                                        BVxz1(:,m+1)=LmVxz1*S(index+m+1);

                                        if m==1
                                            %Vxz
                                            AVxz2(:,m+1)=LmVxz2*C(index+m-1);
                                            BVxz2(:,m+1)=LmVxz2*S(index+m-1);
                                                
                                            %Vyz
                                            AVyz2(:,m+1)=LmVyz2*C(index+m-1).*q;
                                            BVyz2(:,m+1)=LmVyz2*S(index+m-1).*q;
                                        else
                                            %Vxz
                                            AVxz2(:,m+1)=0;
                                            BVxz2(:,m+1)=0;
                                                
                                            %Vyz
                                            AVyz2(:,m+1)=0;
                                            BVyz2(:,m+1)=0;
                                        end
                                        
                                        %Vxy
                                        AVxy1(:,m+1)=(LmVxy1*C(index+m+2)).*q;  
                                        BVxy1(:,m+1)=(LmVxy1*S(index+m+2)).*q; 
                                           
                                        AVxy3(:,m+1)=LmVxy3*zeros(nmax-2+1,1);
                                        BVxy3(:,m+1)=LmVxy3*zeros(nmax-2+1,1);
                                            
                                        %Vyz
                                        AVyz1(:,m+1)=(LmVyz1*C(index+m+1)).*q;  
                                        BVyz1(:,m+1)=(LmVyz1*S(index+m+1)).*q; 
                                    elseif m>nmax-2
                                            
                                        if m==nmax
                                            %Vxz
                                            AVxz1(:,m+1)=0;
                                            BVxz1(:,m+1)=0;
                                                
                                            %Vyz
                                            AVyz1(:,m+1)=0;
                                            BVyz1(:,m+1)=0;
                                        else
                                            %Vxz
                                            AVxz1(:,m+1)=LmVxz1*C(index((m-1):end)+m+1);
                                            BVxz1(:,m+1)=LmVxz1*S(index((m-1):end)+m+1);
                                                
                                            %Vyz
                                            AVyz1(:,m+1)=LmVyz1*C(index((m-1):end)+m+1).*q;
                                            BVyz1(:,m+1)=LmVyz1*S(index((m-1):end)+m+1).*q;
                                        end
                                           
                                        AVxz2(:,m+1)=LmVxz2*C(index((m-1):end)+m-1);
                                        BVxz2(:,m+1)=LmVxz2*S(index((m-1):end)+m-1);

                                        %Vxy
                                        AVxy1(:,m+1)=0;
                                        BVxy1(:,m+1)=0;
                                            
                                        AVxy3(:,m+1)=(LmVxy3*C(index((m-1):end)+m-2)).*q;
                                        BVxy3(:,m+1)=(LmVxy3*S(index((m-1):end)+m-2)).*q;  
                                            
                                        %Vyz
                                        AVyz2(:,m+1)=LmVyz2*C(index((m-1):end)+m-1).*q;
                                        BVyz2(:,m+1)=LmVyz2*S(index((m-1):end)+m-1).*q;
                                    else
                                        AVxz1(:,m+1)=LmVxz1*C(index((m-1):end)+m+1);
                                        BVxz1(:,m+1)=LmVxz1*S(index((m-1):end)+m+1);

                                        AVxz2(:,m+1)=LmVxz2*C(index((m-1):end)+m-1);
                                        BVxz2(:,m+1)=LmVxz2*S(index((m-1):end)+m-1);

                                        %Vxy
                                        AVxy1(:,m+1)=(LmVxy1*C(index((m-1):end)+m+2)).*q;
                                        BVxy1(:,m+1)=(LmVxy1*S(index((m-1):end)+m+2)).*q;
                                            
                                        AVxy3(:,m+1)=(LmVxy3*C(index((m-1):end)+m-2)).*q;
                                        BVxy3(:,m+1)=(LmVxy3*S(index((m-1):end)+m-2)).*q;
                                            
                                        %Vyz
                                        AVyz1(:,m+1)=LmVyz1*C(index((m-1):end)+m+1).*q;
                                        BVyz1(:,m+1)=LmVyz1*S(index((m-1):end)+m+1).*q;
                                           
                                        AVyz2(:,m+1)=LmVyz2*C(index((m-1):end)+m-1).*q;
                                        BVyz2(:,m+1)=LmVyz2*S(index((m-1):end)+m-1).*q;
                                    end
                                elseif radenie==1
                                    if m<2
                                        %Vxy
                                        AVxy1(:,m+1)=(LmVxy1*[zeros(m,1);C((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))]).*q;
                                        BVxy1(:,m+1)=(LmVxy1*[zeros(m,1);S((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))]).*q;

                                        AVxy3(:,m+1)=0;
                                        BVxy3(:,m+1)=0;
                                            
                                        %Vxz
                                        AVxz1(:,m+1)=(LmVxz1*C((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)))); 
                                        BVxz1(:,m+1)=(LmVxz1*S((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)))); 
                                            
                                        %Vyz
                                        AVyz1(:,m+1)=(LmVyz1*C((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)))).*q; 
                                        BVyz1(:,m+1)=(LmVyz1*S((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)))).*q; 
                                            
                                        if m==1
                                            %Vxz
                                            AVxz2(:,m+1)=LmVxz2*C(1:(nmax-1));
                                            BVxz2(:,m+1)=0;
                                                
                                            %Vyz
                                            AVyz2(:,m+1)=(LmVyz2*C(1:(nmax-1))).*q;
                                            BVyz2(:,m+1)=0;
                                        else
                                            %Vxz
                                            AVxz2(:,m+1)=0;
                                            BVxz2(:,m+1)=0;
                                                
                                            %Vyz
                                            AVyz2(:,m+1)=0;
                                            BVyz2(:,m+1)=0;
                                        end
                                    elseif m>nmax-2
                                        if m==2
                                            %Vxy
                                            AVxy3(:,m+1)=(LmVxy3*C(1:(nmax-1))).*q;
                                            BVxy3(:,m+1)=0;
                                        else
                                            %Vxy
                                            AVxy3(:,m+1)=(LmVxy3*C((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)))).*q;
                                            BVxy3(:,m+1)=(LmVxy3*S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)))).*q; 
                                        end   
                                            
                                        AVxy1(:,m+1)=0;
                                        BVxy1(:,m+1)=0;
                                            
                                        if m==nmax
                                            %Vxz
                                            AVxz1(:,m+1)=0;    
                                            BVxz1(:,m+1)=0;
                                                
                                            %Vyz
                                            AVyz1(:,m+1)=0;    
                                            BVyz1(:,m+1)=0;
                                        else
                                            %Vxz
                                            AVxz1(:,m+1)=LmVxz1*[0;C((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                            BVxz1(:,m+1)=LmVxz1*[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                                
                                            %Vyz
                                            AVyz1(:,m+1)=(LmVyz1*[0;C((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]).*q;
                                            BVyz1(:,m+1)=(LmVyz1*[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]).*q;
                                        end

                                        AVxz2(:,m+1)=LmVxz2*C((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                        BVxz2(:,m+1)=LmVxz2*S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                          
                                        %Vyz
                                        AVyz2(:,m+1)=LmVyz2*C((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))).*q;
                                        BVyz2(:,m+1)=LmVyz2*S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))).*q;
                                    else                   
                                        %Vxy
                                        AVxy1(:,m+1)=(LmVxy1*[0;0;C((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))]).*q;
                                        BVxy1(:,m+1)=(LmVxy1*[0;0;S((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))]).*q;

                                        if m==2
                                            AVxy3(:,m+1)=(LmVxy3*C(1:(nmax-1))).*q;
                                            BVxy3(:,m+1)=0;
                                        else
                                            AVxy3(:,m+1)=(LmVxy3*C((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)))).*q;
                                            BVxy3(:,m+1)=(LmVxy3*S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)))).*q;
                                        end
                                        
                                        %Vxz
                                        AVxz1(:,m+1)=(LmVxz1*[0;C((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]);
                                        BVxz1(:,m+1)=(LmVxz1*[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]);
                                            
                                        AVxz2(:,m+1)=LmVxz2*C((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                        BVxz2(:,m+1)=LmVxz2*S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                            
                                        %Vyz
                                        AVyz2(:,m+1)=LmVyz2*C((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))).*q;
                                        BVyz2(:,m+1)=LmVyz2*S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))).*q;
                                            
                                        AVyz1(:,m+1)=(LmVyz1*[0;C((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]).*q;
                                        BVyz1(:,m+1)=(LmVyz1*[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))]).*q;
                                    end
                                end
                                
                                %Modified forward column method cannot be
                                %applied
                                
                            elseif volbapar(i)==16 %Gravity  

                                if m==nmax                                
                                    ampl_Wr=((2:nmax)+1);
                                end

                                if m<2
                                    LmWr=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_Wr);
                                    LmWlambda=m*Pnm(:,(3-m):(end-m));
                                    LmWfi=dPnm(:,(3-m):(end-m));
                                else
                                    LmWr=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Wr(:,(m-1):end));
                                    LmWlambda=m*Pnm(:,1:(nmax-m+1));
                                    LmWfi=dPnm(:,1:(nmax-m+1));
                                end                           

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AWr=zeros(length(fiG),nmax+1);
                                        BWr=zeros(length(fiG),nmax+1);
                                        AWlambda=AWr;
                                        BWlambda=BWr;
                                        AWfi=AWr;
                                        BWfi=BWr;
                                    end

                                    AWr(:,m+1)=LmWr*Cm;
                                    BWr(:,m+1)=LmWr*Sm;   
                                    AWlambda(:,m+1)=LmWlambda*Cm;
                                    BWlambda(:,m+1)=LmWlambda*Sm;
                                    AWfi(:,m+1)=LmWfi*Cm;
                                    BWfi(:,m+1)=LmWfi*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Wr=zeros(length(fiG),length(lambda));
                                        Wlambda=Wr;
                                        Wfi=Wr;
                                    end

                                    Wr=bsxfun(@times,Wr,u)+(LmWr*Cm*cos(m*lambda')+LmWr*Sm*sin(m*lambda'));
                                    Wlambda=bsxfun(@times,Wlambda,u)+(-LmWlambda*Cm*sin(m*lambda')+LmWlambda*Sm*cos(m*lambda'));
                                    Wfi=bsxfun(@times,Wfi,u)+(LmWfi*Cm*cos(m*lambda')+LmWfi*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==17 %Gravity sa

                                if m==nmax                                
                                    ampl_g_sa=((2:nmax)+1);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_g_sa);
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_g_sa(:,(m-1):end));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        Ag_sa=zeros(length(fiG),nmax+1);
                                        Bg_sa=zeros(length(fiG),nmax+1);
                                    end

                                    Ag_sa(:,m+1)=Lm*Cm;
                                    Bg_sa(:,m+1)=Lm*Sm;                             
                                elseif volbaALFs==2
                                    if m==nmax
                                       g_sa=zeros(length(fiG),length(lambda)); 
                                    end

                                    g_sa=bsxfun(@times,g_sa,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==18 %Gravity potential                           

                                if m<2
                                    Lm=Pnm(:,(3-m):(end-m));
                                else
                                    Lm=Pnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3                               
                                    if m==nmax
                                        AW=zeros(length(fiG),nmax+1);
                                        BW=zeros(length(fiG),nmax+1);
                                    end

                                    AW(:,m+1)=Lm*Cm;
                                    BW(:,m+1)=Lm*Sm;                    
                                elseif volbaALFs==2
                                    if m==nmax
                                        W=zeros(length(fiG),length(lambda));
                                    end

                                    W=bsxfun(@times,W,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==19 %Gravity anomaly sa

                                if m==nmax
                                    ampl_anomalia_sa=((2:nmax)-1);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_anomalia_sa);
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_anomalia_sa(:,(m-1):end));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        Aanomalia_sa=zeros(length(fiG),nmax+1);
                                        Banomalia_sa=zeros(length(fiG),nmax+1);
                                    end

                                    Aanomalia_sa(:,m+1)=Lm*deltaCm;
                                    Banomalia_sa(:,m+1)=Lm*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        anomalia_sa=zeros(length(fiG),length(lambda));
                                    end

                                    anomalia_sa=bsxfun(@times,anomalia_sa,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==20 %Gravity disturbance

                                if m==nmax                                
                                    %ampl_Wr=zeros(length(fiG),nmax-1); %Damping factor
                                    %for n=2:nmax
                                    %    ampl_Wr(:,n-1)=(n+1);
                                    %end
                                    ampl_Wrpor=((2:nmax)+1);
                                end

                                if m<2
                                    %LmWr=Pnm(:,(3-m):(end-m)).*ampl_Wr;
                                    LmWrpor=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_Wrpor);
                                    LmWlambdapor=m*Pnm(:,(3-m):(end-m));
                                    LmWfipor=dPnm(:,(3-m):(end-m));
                                else
                                    %LmWr=Pnm(:,1:(nmax-m+1)).*ampl_Wr(:,(m-1):end);
                                    LmWrpor=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Wrpor(:,(m-1):end));
                                    LmWlambdapor=m*Pnm(:,1:(nmax-m+1));
                                    LmWfipor=dPnm(:,1:(nmax-m+1));
                                end                           

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AWrpor=zeros(length(fiG),nmax+1);
                                        BWrpor=zeros(length(fiG),nmax+1);
                                        AWlambdapor=AWrpor;
                                        BWlambdapor=BWrpor;
                                        AWfipor=AWrpor;
                                        BWfipor=BWrpor;
                                    end

                                    AWrpor(:,m+1)=LmWrpor*Cm;
                                    BWrpor(:,m+1)=LmWrpor*Sm;   
                                    AWlambdapor(:,m+1)=LmWlambdapor*Cm;
                                    BWlambdapor(:,m+1)=LmWlambdapor*Sm;
                                    AWfipor(:,m+1)=LmWfipor*Cm;
                                    BWfipor(:,m+1)=LmWfipor*Sm;

                                    if m==0
                                        AUr(:,m+1)=LmWrpor*CElm;
                                        AUfi(:,m+1)=LmWfipor*CElm;
                                    end
                                elseif volbaALFs==2
                                    if m==nmax
                                        Wrpor=zeros(length(fiG),length(lambda));
                                        Wlambdapor=Wrpor;
                                        Wfipor=Wrpor;
                                        Ur=Wrpor;
                                        Ufi=Ur;
                                    end

                                    Wrpor=bsxfun(@times,Wrpor,u)+(LmWrpor*Cm*cos(m*lambda')+LmWrpor*Sm*sin(m*lambda'));
                                    Wlambdapor=bsxfun(@times,Wlambdapor,u)+(-LmWlambdapor*Cm*sin(m*lambda')+LmWlambdapor*Sm*cos(m*lambda'));
                                    Wfipor=bsxfun(@times,Wfipor,u)+(LmWfipor*Cm*cos(m*lambda')+LmWfipor*Sm*sin(m*lambda'));

                                    if m==0                                   
                                        Ur=LmWrpor*CElm*cos(m*lambda');
                                        Ufi=LmWfipor*CElm*cos(m*lambda');
                                    end
                                end

                            elseif volbapar(i)==21 %Gravity disturbance sa

                                if m==nmax                                
                                    ampl_porucha_sa=((2:nmax)+1);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_porucha_sa);
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_porucha_sa(:,(m-1):end));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        Aporucha_sa=zeros(length(fiG),nmax+1);
                                        Bporucha_sa=zeros(length(fiG),nmax+1);
                                    end

                                    Aporucha_sa(:,m+1)=Lm*deltaCm;
                                    Bporucha_sa(:,m+1)=Lm*Sm;                            
                                elseif volbaALFs==2
                                    if m==nmax
                                        porucha_sa=zeros(length(fiG),length(lambda));
                                    end

                                    porucha_sa=bsxfun(@times,porucha_sa,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==22 %Height anomaly ell

                                if m<2
                                    Lm=Pnm(:,(3-m):(end-m));
                                else
                                    Lm=Pnm(:,1:(nmax-m+1));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AzetaEl=zeros(length(fiG),nmax+1);
                                        BzetaEl=zeros(length(fiG),nmax+1);
                                    end

                                    AzetaEl(:,m+1)=Lm*deltaCm;
                                    BzetaEl(:,m+1)=Lm*Sm;                            
                                elseif volbaALFs==2
                                    if m==nmax
                                        zetaEl=zeros(length(fiG),length(lambda));
                                    end

                                    zetaEl=bsxfun(@times,zetaEl,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==23 %Height anomaly
                                
                                if m==nmax
                                    ampl_zeta_H=zeros(length(fiG),nmax+1); %Damping factor
                                    for n=0:nmax
                                        ampl_zeta_H(:,n+1)=1./((R./r).^n); 
                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFFIED
                                        % fnALFS
                                    end
                                    
                                    ampl_zeta_dg=((2:nmax)+1);
                                end
                                        
                                if m<2
                                    Lmdg=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_zeta_dg);
                                    Lm=Pnm(:,(3-m):(end-m));
                                    LmH=Pnm(:,1:(nmax-m+1)).*ampl_zeta_H(:,(m+1):end);
                                else
                                    Lmdg=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_zeta_dg(:,(m-1):end));
                                    Lm=Pnm(:,1:(nmax-m+1));
                                    LmH=Pnm(:,1:(nmax-m+1)).*ampl_zeta_H(:,(m+1):end);
                                end
                                
                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AN1c_zeta=zeros(length(fiG),nmax);
                                        BN1c_zeta=AN1c_zeta;
                                        Azetadg=AN1c_zeta;
                                        Bzetadg=AN1c_zeta;
                                        AH_zeta=AN1c_zeta;
                                        BH_zeta=AN1c_zeta;
                                        Azeta=AN1c_zeta;
                                        Bzeta=AN1c_zeta;
                                    end  
                                    
                                    AN1c_zeta(:,m+1)=Lm*deltaCm;
                                    BN1c_zeta(:,m+1)=Lm*Sm;  
                                    Azetadg(:,m+1)=Lmdg*deltaCm;
                                    Bzetadg(:,m+1)=Lmdg*Sm;
                                    AH_zeta(:,m+1)=LmH*HCm;
                                    BH_zeta(:,m+1)=LmH*HSm;
                                    Azeta(:,m+1)=Lm*deltaCm;
                                    Bzeta(:,m+1)=Lm*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        zeta_N1c=zeros(length(fiG),length(lambda));
                                        zeta_H=zeta_N1c;
                                        zeta_dg=zeta_N1c;
                                        zeta_zetaEl=zeta_N1c;
                                    end

                                    zeta_N1c=bsxfun(@times,zeta_N1c,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                    zeta_H=bsxfun(@times,zeta_H,u)+(LmH*HCm*cos(m*lambda')+LmH*HSm*sin(m*lambda'));
                                    zeta_dg=bsxfun(@times,zeta_dg,u)+(Lmdg*deltaCm*cos(m*lambda')+Lmdg*Sm*sin(m*lambda'));
                                    zeta_zetaEl=bsxfun(@times,zeta_zetaEl,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end
                                
                                if m==0
                                    clear ampl_zeta_H ampl_zeta_dg
                                end
                                                                
                            elseif volbapar(i)==24 %Second radial derivative of disturbing potential

                                if m==nmax                               
                                    ampl_T_rr=((2:nmax)+1).*((2:nmax)+2);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_T_rr);
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_T_rr(:,(m-1):end));
                                end
                                
                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AT_rr=zeros(length(fiG),nmax+1);
                                        BT_rr=zeros(length(fiG),nmax+1);
                                    end

                                    AT_rr(:,m+1)=Lm*deltaCm;
                                    BT_rr(:,m+1)=Lm*Sm;                             
                                elseif volbaALFs==2
                                    if m==nmax
                                        T_rr=zeros(length(fiG),length(lambda));
                                    end

                                    T_rr=bsxfun(@times,T_rr,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==25 %Second radial derivative of disturbing potential

                                if m==nmax                               
                                    ampl_Wrr=((2:nmax)+1).*((2:nmax)+2);
                                end

                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,(3-m):(end-m)),ampl_Wrr);
                                else
                                    Lm=bsxfun(@times,Pnm(:,1:(nmax-m+1)),ampl_Wrr(:,(m-1):end));
                                end

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AWrr=zeros(length(fiG),nmax+1);
                                        BWrr=zeros(length(fiG),nmax+1);
                                    end

                                    AWrr(:,m+1)=Lm*Cm;
                                    BWrr(:,m+1)=Lm*Sm;                           
                                elseif volbaALFs==2
                                    if m==nmax
                                        Wrr=zeros(length(fiG),length(lambda));
                                    end

                                    Wrr=bsxfun(@times,Wrr,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end
                            end
                        end                                       
                    end

                    %Update of the progress bar
                    set(progressbar,'string',...
                        'Progress: Matrix multiplications...',...
                        'fontsize',8); drawnow;

                    clear Lm dLm Pnm dPnm ddPnm Cm Sm C CEl CElm deltaC ...
                        deltaCm S u t q q2 index tu qu enm singdPnm singddPnm
                    
                    clear ampl_Trr ampl_Tfl ampl_Tzz amplH ampl_Vrr ampl_Vfl ...
                        ampl_Wr ampl_g_sa ampl_anomalia_sa ampl_Wrpor ...
                        ampl_porucha_sa ampl_T_rr ampl_Wrr
                    
                    if volbaALFs==1 || volbaALFs==3
                        cosla=cos((0:nmax)'*lambda');
                        sinla=sin((0:nmax)'*lambda');
                    end

                    %Computation of the normal gravity for eta, xi, Theta,
                    %Geoid undulation, zeta el, zeta
                    if any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==10) || any(volbapar==22) || any(volbapar==23)                        
                        bEl=aEl*sqrt(1-eEl^2);
                        EEl=sqrt(aEl^2-bEl^2);

                        %Computation of ellipsoidal harmonic coordinates
                        [X,Y,Z]=geodetic2ecef(fi,0*zeros(length(fi),1),h,[aEl eEl]');

                        ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                        betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                        clear X Y Z
                        
                        wgama=sqrt((ugama.^2+EEl^2*sin(betagama).^2)./(ugama.^2+EEl^2));
                        qgama=1/2*((1+(3*ugama.^2)./EEl^2).*atan(EEl./ugama)-3*ugama./EEl);
                        qgama_=3*(1+(ugama.^2)./EEl^2).*(1-ugama./EEl.*atan(EEl./ugama))-1;
                        qgama0=1/2*((1+(3*bEl^2)/EEl^2)*atan(EEl/bEl)-3*bEl/EEl);
                        
                        gamau=-1./wgama.*(GMEl./(ugama.^2+EEl.^2)+(omegaEl.^2.*aEl.^2.*EEl)./(ugama.^2+EEl.^2).*(qgama_./qgama0).*(1./2.*sin(betagama).^2-1./6)-omegaEl.^2.*ugama.*cos(betagama).^2);
                        gamabeta=-1./wgama.*(-(omegaEl.^2.*aEl.^2.*qgama)./(sqrt(ugama.^2+EEl.^2).*qgama0)+omegaEl.^2.*sqrt(ugama.^2+EEl.^2)).*sin(betagama).*cos(betagama);
                        
                        clear ugama betagama wgama qgama qgama_ qgama0
                        
                        gamaP=sqrt(gamau.^2+gamabeta.^2);
                        
                        clear gamau gamabeta
                    end

                    %% Final computation of functionals of the geopotential             
                    for i=1:pocetpar
                        if volbapar(i)==1                
                        elseif volbapar(i)==2 %Deflection of the vertical eta
                            if volbaALFs==1 || volbaALFs==3
                                eta=-Aeta*sinla+Beta*cosla;                                
                                clear Aeta Beta
                            elseif volbaALFs==2
                                eta=eta*1e280;
                            end

                            Pg=bsxfun(@times,-GM./(r.^2.*gamaP.*cos(fiG)),(eta))*(180/pi)*3600;
                            Pg(fi==pi/2 | fi==-pi/2,:)=0;
                            Pg=Pg(:);

                            clear eta
                        elseif volbapar(i)==3 %Deflection of the vertical xi
                            if volbaALFs==1 || volbaALFs==3
                                ksi=Aksi*cosla+Bksi*sinla;
                                clear Aksi Bksi
                            elseif volbaALFs==2
                                ksi=ksi*1e280;
                            end

                            Pg=bsxfun(@times,-GM./(r.^2.*gamaP),(ksi))*(180/pi)*3600;
                            Pg=Pg(:);

                            clear ksi
                        elseif volbapar(i)==4 %Deflection of the vertical Theta
                            if volbaALFs==1 || volbaALFs==3
                                Teta=-ATeta*sinla+BTeta*cosla;
                                clear ATeta BTeta
                                Tksi=ATksi*cosla+BTksi*sinla;
                                clear ATksi BTksi
                            elseif volbaALFs==2
                                Teta=Teta*1e280;
                                Tksi=Tksi*1e280;
                            end

                            Teta=bsxfun(@times,-GM./(r.^2.*gamaP.*cos(fiG)),(Teta))*(180/pi)*3600;
                            Teta(fi==pi/2 | fi==-pi/2,:)=0;
                            Tksi=bsxfun(@times,-GM./(r.^2.*gamaP),(Tksi))*(180/pi)*3600;
                            Teta=Teta(:);
                            Tksi=Tksi(:);
                            Talfa=atan2(Teta,Tksi);
                            Talfa(Talfa<0)=Talfa(Talfa<0)+2*pi;
   
                            Pg=[sqrt(Teta.^2+Tksi.^2) rad2deg(Talfa)];

                            clear Teta Tksi Talfa
                        elseif volbapar(i)==5 %Disturbing potential                        
                            if volbaALFs==1 || volbaALFs==3
                                T=AT*cosla+BT*sinla;
                                clear AT BT
                            elseif volbaALFs==2
                                T=T*1e280;
                            end

                            Pg=bsxfun(@times,GM./r,T);
                            Pg=Pg(:);

                            clear T
                        elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                            
                            clear Lmff Lmll                            
                            
                            if volbaALFs==1 || volbaALFs==3
                                Trr=ATrr*cosla+BTrr*sinla;
                                clear ATrr BTrr
                                Tff=ATff*cosla+BTff*sinla;
                                clear ATff BTff
                                
                                ATll(fi>deg2rad(89.9) | fi<deg2rad(-89.9),:)=0;
                                BTll(fi>deg2rad(89.9) | fi<deg2rad(-89.9),:)=0;
                                Tll=ATll*cosla+BTll*sinla;
                                clear ATll BTll
                            elseif volbaALFs==2
                                Trr=Trr*1e280;
                                Tff=Tff*1e280;
                                
                                Tll(fi>deg2rad(89.9) | fi<deg2rad(-89.9),:)=0;
                                Tll=Tll*1e280;
                            end

                            Trr=bsxfun(@times,GM./r.^3,Trr)*10^9;
                            Tff=bsxfun(@times,GM./r.^3,Tff)*10^9;
                            Tll=bsxfun(@times,GM./(r.^3.*cos(fiG).^2),Tll)*10^9;
                            
                            Pg=Trr(:);
                            clear Trr
                            Pg=[Pg Tff(:)];
                            clear Tff
                            Pg=[Pg -Tll(:)];
                            clear Tll  
                            
                        elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                            
                            clear Lmrf Lmrl Lmfl
                            
                            if volbaALFs==1 || volbaALFs==3
                                Trf=ATrf*cosla+BTrf*sinla;
                                clear ATrf BTrf
                                Trl=-ATrl*sinla+BTrl*cosla;
                                clear ATrl BTrl
                                
                                ATfl(fi>deg2rad(89.5) | fi<deg2rad(-89.5),:)=0;
                                BTfl(fi>deg2rad(89.5) | fi<deg2rad(-89.5),:)=0;
                                Tfl=-ATfl*sinla+BTfl*cosla;
                                clear ATfl BTfl
                            elseif volbaALFs==2
                                Trf=Trf*1e280;
                                Trl=Trl*1e280;
                                
                                Tfl(fi>deg2rad(89.5) | fi<deg2rad(-89.5),:)=0;
                                Tfl=Tfl*1e280;
                            end

                            Trf=bsxfun(@times,GM./r.^3,Trf)*10^9;
                            Trl=bsxfun(@times,GM./(r.^3.*cos(fiG)),Trl)*10^9;
                            Tfl=bsxfun(@times,GM./(r.^3.*cos(fiG)),Tfl)*10^9;
                            
                            Pg=-Trf(:);
                            clear Trf
                            Pg=[Pg -Trl(:)];
                            clear Trl
                            Pg=[Pg Tfl(:)];
                            clear Tfl  
                            
                        elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                            
                            clear anm bnm cnm LmTxx1 LmTxx2 LmTxx3 ...
                                LmTyy1 ampl_Tzz
                            
                            if volbaALFs==1 || volbaALFs==3
                                Tzz=ATzz*cosla+BTzz*sinla;
                                clear ATzz BTzz                                
                                
                                %Txx
                                Txx1=ATxx1*cosla+BTxx1*sinla;
                                clear ATxx1 BTxx1
                                
                                Txx2=ATxx2*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]+BTxx2*[sinla(3:end,:);zeros(2,length(cosla(1,:)))];
                                clear ATxx2 BTxx2
                                
                                Txx=Txx1+Txx2;
                                
                                %Tyy
                                Tyy1=ATyy1*cosla+BTyy1*sinla;
                                clear ATyy1 BTyy1
                                
                                Tyy=Tyy1+Txx2;
                                clear Txx1 Txx2 Tyy1
                                                         
                                Txx3=ATxx3*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]+BTxx3*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)];
                                clear ATxx3 BTxx3
                                
                                Txx=Txx+Txx3;                                                                                                                       
                                Tyy=Tyy+Txx3;
                                clear Txx3
                                
                            elseif volbaALFs==2
                            end

                            Tzz=bsxfun(@times,GM./r.^3,Tzz)*10^9;
                            
                            Txx=bsxfun(@times,GM./r.^3,Txx)*10^9;
                            
                            Tyy=bsxfun(@times,GM./r.^3,Tyy)*10^9;
                            
                            Pg=Txx(:);
                            clear Txx
                            Pg=[Pg -Tyy(:)];
                            clear Tyy
                            Pg=[Pg Tzz(:)];
                            clear Tzz

                        elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                            
                            clear gamanm betanm gnm hnm dnm minm ninm ...
                                LmTxz1 LmTxz2 LmTxy1 LmTxy2 LmTxy3 ...
                                LmTyz1 LmTyz2
                            
                            if volbaALFs==1 || volbaALFs==3

                                %Txz
                                Txz1=ATxz1*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]+BTxz1*[sinla(2:end,:);zeros(1,length(cosla(1,:)))];
                                clear ATxz1 BTxz1
                                
                                Txz2=ATxz2*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]+BTxz2*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)];
                                clear ATxz2 BTxz2
                                
                                Txz=Txz1+Txz2;
                                clear Txz1 Txz2
                                
                                %Txy
                                Txy1=-ATxy1*[sinla(3:end,:);zeros(2,length(sinla(1,:)))]+BTxy1*[cosla(3:end,:);zeros(2,length(cosla(1,:)))];
                                clear ATxy1 BTxy1
                                
                                Txy2=-ATxy2*sinla+BTxy2*cosla;
                                clear ATxy2 BTxy2
                                
                                Txy3=-ATxy3*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]+BTxy3*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)];
                                clear ATxy3 BTxy3

                                Txy=Txy1+Txy2+Txy3;
                                clear Txy1 Txy2 Txy3
                                
                                %Tyz
                                Tyz1=-ATyz1*[sinla(2:end,:);zeros(1,length(sinla(1,:)))]+BTyz1*[cosla(2:end,:);zeros(1,length(cosla(1,:)))];
                                clear ATyz1 BTyz1
                                
                                Tyz2=-ATyz2*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]+BTyz2*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)];
                                clear ATyz2 BTyz2
 
                                Tyz=Tyz1+Tyz2;
                                clear Tyz1 Tyz2
                            elseif volbaALFs==2
                            end

                            Txz=bsxfun(@times,GM./r.^3,Txz)*10^9;
                            Txy=bsxfun(@times,GM./r.^3,Txy)*10^9;
                            Tyz=bsxfun(@times,GM./r.^3,Tyz)*10^9;

                            Pg=Txy(:);
                            clear Txy
                            Pg=[Pg Txz(:)];
                            clear Txz
                            Pg=[Pg Tyz(:)];
                            clear Tyz

                        elseif volbapar(i)==10 %Geoid undulation
                            
                            clear HC HCm HS HSm LmH
                            
                            if volbaALFs==1 || volbaALFs==3
                                H=AH*cosla+BH*sinla;
                                clear AH BH
                                N1c=AN1c*cosla+BN1c*sinla;
                                clear AN1c BN1c
                            elseif volbaALFs==2
                                H=H*1e280;
                                N1c=N1c*1e280;
                            end                            
                            
                            N1c=bsxfun(@times,GM./(r.*gamaP),N1c);
                            H(H<0)=H(H<0)*0; %H is set to zero in the areas of oceans and seas
                           
                            G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                            ro=2670; %Density of the crust

                            Pg=N1c-bsxfun(@times,(2*pi*G*ro*H.^2),1./gamaP);
                            Pg=Pg(:);
                            
                            clear H N1c
                        elseif volbapar(i)==11 %Gravitational potential
                            if volbaALFs==1 || volbaALFs==3
                                V=AV*cosla+BV*sinla;
                                clear AV BV
                            elseif volbaALFs==2
                                V=V*1e280;
                            end

                            Pg=bsxfun(@times,GM./r,V);

                            if nmin0==1 %Zero degree term
                                Pg=bsxfun(@plus,GM./r,Pg);
                            end
                            
                            Pg=Pg(:);

                            clear V 
                        elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                            
                            clear Lmff Lmll                            
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vrr=AVrr*cosla+BVrr*sinla;
                                clear AVrr BVrr
                                Vff=AVff*cosla+BVff*sinla;
                                clear AVff BVff
                                
                                AVll(fi>deg2rad(89.9) | fi<deg2rad(-89.9),:)=0;
                                BVll(fi>deg2rad(89.9) | fi<deg2rad(-89.9),:)=0;
                                Vll=AVll*cosla+BVll*sinla;
                                clear AVll BVll
                            elseif volbaALFs==2
                                Vrr=Vrr*1e280;
                                Vff=Vff*1e280;
                                
                                Vll(fi>deg2rad(89.9) | fi<deg2rad(-89.9),:)=0;
                                Vll=Vll*1e280;
                            end

                            Vrr=bsxfun(@times,GM./r.^3,Vrr);
                            if nmin0==1 %Zero degree term
                                Vrr=bsxfun(@plus,2*GM./r.^3,Vrr);
                            end
                            Vrr=Vrr*10^9;
                            
                            Vff=bsxfun(@times,GM./r.^3,Vff)*10^9;
                            Vll=bsxfun(@times,GM./(r.^3.*cos(fiG).^2),Vll)*10^9;
                            
                            Pg=Vrr(:);
                            clear Vrr
                            Pg=[Pg Vff(:)];
                            clear Vff
                            Pg=[Pg -Vll(:)];
                            clear Vll  
                            
                        elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                            
                            clear Lmrf Lmrl Lmfl
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vrf=AVrf*cosla+BVrf*sinla;
                                clear AVrf BVrf
                                Vrl=-AVrl*sinla+BVrl*cosla;
                                clear AVrl BVrl
                                
                                AVfl(fi>deg2rad(89.5) | fi<deg2rad(-89.5),:)=0;
                                BVfl(fi>deg2rad(89.5) | fi<deg2rad(-89.5),:)=0;
                                Vfl=-AVfl*sinla+BVfl*cosla;
                                clear AVfl BVfl
                            elseif volbaALFs==2
                                Vrf=Vrf*1e280;
                                Vrl=Vrl*1e280;
                                
                                Vfl(fi>deg2rad(89.5) | fi<deg2rad(-89.5),:)=0;
                                Vfl=Vfl*1e280;
                            end

                            Vrf=bsxfun(@times,GM./r.^3,Vrf)*10^9;
                            Vrl=bsxfun(@times,GM./(r.^3.*cos(fiG)),Vrl)*10^9;
                            Vfl=bsxfun(@times,GM./(r.^3.*cos(fiG)),Vfl)*10^9;
                            
                            Pg=-Vrf(:);
                            clear Vrf
                            Pg=[Pg -Vrl(:)];
                            clear Vrl
                            Pg=[Pg Vfl(:)];
                            clear Vfl  
                            
                        elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                            
                            clear anm bnm cnm LmVxx1 LmVxx2 LmVxx3 ...
                                LmVyy1 ampl_Vzz
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vzz=AVzz*cosla+BVzz*sinla;
                                clear AVzz BVzz                                
                                
                                %Vxx
                                Vxx1=AVxx1*cosla+BVxx1*sinla;
                                clear AVxx1 BVxx1
                                
                                Vxx2=AVxx2*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]+BVxx2*[sinla(3:end,:);zeros(2,length(cosla(1,:)))];
                                clear AVxx2 BVxx2
                                
                                Vxx=Vxx1+Vxx2;
                                
                                %Vyy
                                Vyy1=AVyy1*cosla+BVyy1*sinla;
                                clear AVyy1 BVyy1
                                
                                Vyy=Vyy1+Vxx2;
                                clear Vxx1 Vxx2 Vyy1
                                                         
                                Vxx3=AVxx3*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]+BVxx3*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)];
                                clear AVxx3 BVxx3
                                
                                Vxx=Vxx+Vxx3;                                                                                                                       
                                Vyy=Vyy+Vxx3;
                                clear Vxx3
                                
                            elseif volbaALFs==2
                            end

                            if nmin0==1 %Zero degree term
                                Vzz=bsxfun(@times,GM./r.^3,2+Vzz)*10^9;
                            elseif nmin0==0
                                Vzz=bsxfun(@times,GM./r.^3,Vzz)*10^9;
                            end
                            
                            if nmin0==1 %Zero degree term
                                Vxx=bsxfun(@times,GM./r.^3,-1+Vxx)*10^9;
                            elseif nmin0==0
                                Vxx=bsxfun(@times,GM./r.^3,Vxx)*10^9;
                            end
                            
                            if nmin0==1
                                Vyy=bsxfun(@times,GM./r.^3,1+Vyy)*10^9;
                            elseif nmin0==0
                                Vyy=bsxfun(@times,GM./r.^3,Vyy)*10^9;
                            end
                            
                            Pg=Vxx(:);
                            clear Vxx
                            Pg=[Pg -Vyy(:)];
                            clear Vyy
                            Pg=[Pg Vzz(:)];
                            clear Vzz
                            
                        elseif volbapar(i)==15 % Gravitational tensor Vxy_Vxz_Vyz
                            
                            clear gamanm betanm gnm hnm dnm minm ninm ...
                                LmVxz1 LmVxz2 LmVxy1 LmVxy2 LmVxy3 ...
                                LmVyz1 LmVyz2
                            
                            if volbaALFs==1 || volbaALFs==3

                                %Vxz
                                Vxz1=AVxz1*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]+BVxz1*[sinla(2:end,:);zeros(1,length(cosla(1,:)))];
                                clear AVxz1 BVxz1
                                
                                Vxz2=AVxz2*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]+BVxz2*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)];
                                clear AVxz2 BVxz2
                                
                                Vxz=Vxz1+Vxz2;
                                clear Vxz1 Vxz2
                                
                                %Vxy
                                Vxy1=-AVxy1*[sinla(3:end,:);zeros(2,length(sinla(1,:)))]+BVxy1*[cosla(3:end,:);zeros(2,length(cosla(1,:)))];
                                clear AVxy1 BVxy1
                                
                                Vxy2=-AVxy2*sinla+BVxy2*cosla;
                                clear AVxy2 BVxy2
                                
                                Vxy3=-AVxy3*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]+BVxy3*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)];
                                clear AVxy3 BVxy3

                                Vxy=Vxy1+Vxy2+Vxy3;
                                clear Vxy1 Vxy2 Vxy3
                                
                                %Vyz
                                Vyz1=-AVyz1*[sinla(2:end,:);zeros(1,length(sinla(1,:)))]+BVyz1*[cosla(2:end,:);zeros(1,length(cosla(1,:)))];
                                clear AVyz1 BVyz1
                                
                                Vyz2=-AVyz2*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]+BVyz2*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)];
                                clear AVyz2 BVyz2
 
                                Vyz=Vyz1+Vyz2;
                                clear Vyz1 Vyz2
                            elseif volbaALFs==2
                            end

                            Vxz=bsxfun(@times,GM./r.^3,Vxz)*10^9;
                            Vxy=bsxfun(@times,GM./r.^3,Vxy)*10^9;
                            Vyz=bsxfun(@times,GM./r.^3,Vyz)*10^9;

                            Pg=Vxy(:);
                            clear Vxy
                            Pg=[Pg Vxz(:)];
                            clear Vxz
                            Pg=[Pg Vyz(:)];
                            clear Vyz

                        elseif volbapar(i)==16 %Gravity
                            clear LmWr LmWlambda LmWfi
                            
                            if volbaALFs==1 || volbaALFs==3
                                Wr=AWr*cosla+BWr*sinla;
                                clear AWr BWr
                                Wlambda=-AWlambda*sinla+BWlambda*cosla;
                                clear AWlambda BWlambda 
                                Wfi=AWfi*cosla+BWfi*sinla;
                                clear AWfi BWfi 
                            elseif volbaALFs==2
                                Wr=Wr*1e280;
                                Wlambda=Wlambda*1e280;
                                Wfi=Wfi*1e280;
                            end

                            Wr=bsxfun(@times,GM./r.^2,Wr);
                            if nmin0==1 %Zero degree term
                                Wr=bsxfun(@plus,GM./r.^2,Wr);
                            end
                            Wr=bsxfun(@plus,-Wr,omegaEl^2.*r.*(cos(fiG).^2));
                            Wr=Wr(:);

                            Wlambda=bsxfun(@times,GM./r,Wlambda);
                            Wlambda=bsxfun(@times,Wlambda,1./(r.*cos(fiG)));
                            Wlambda=Wlambda(:);

                            Wfi=bsxfun(@times,GM./r,Wfi);
                            Wfi=bsxfun(@plus,Wfi,-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                            Wfi=bsxfun(@times,Wfi,1./r);
                            Wfi=Wfi(:);

                            Pg=sqrt(Wr.^2+Wlambda.^2+Wfi.^2)*10^5;

                            clear Wr Wlambda Wfi
                        elseif volbapar(i)==17 %Gravity sa
                            if volbaALFs==1 || volbaALFs==3
                                g_sa=Ag_sa*cosla+Bg_sa*sinla;
                                clear Ag_sa Bg_sa
                            elseif volbaALFs==2
                                g_sa=g_sa*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^2,g_sa);
                            if nmin0==1 %Zero degree term
                                Pg=bsxfun(@plus,GM./r.^2,Pg);
                            end
                            
                            Pg=sqrt(bsxfun(@plus,-Pg,omegaEl^2.*r.*(cos(fiG).^2)).^2)*10^5;                    
                            Pg=Pg(:);

                            clear g_sa
                        elseif volbapar(i)==18 %Gravity potential
                            if volbaALFs==1 || volbaALFs==3
                                W=AW*cosla+BW*sinla;
                                clear AW BW
                            elseif volbaALFs==2
                                W=W*1e280;
                            end

                            Pg=bsxfun(@times,GM./r,W);
                            clear W
                            if nmin0==1 %Zero degree term
                                Pg=bsxfun(@plus,GM./r,Pg);
                            end
                            Pg=bsxfun(@plus,Pg,1/2*omegaEl.^2.*r.^2.*cos(fiG).^2);
                            Pg=Pg(:);
                        elseif volbapar(i)==19 %Gravity anomaly sa
                            if volbaALFs==1 || volbaALFs==3
                                anomalia_sa=Aanomalia_sa*cosla+Banomalia_sa*sinla;
                                clear Aanomalia_sa Banomalia_sa
                            elseif volbaALFs==2
                                anomalia_sa=anomalia_sa*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^2,anomalia_sa)*10^5;
                            Pg=Pg(:);

                            clear anomalia_sa
                        elseif volbapar(i)==20 %Gravity disturbance
                            clear LmWrpor LmWlambdapor LmWfipor
                            
                            if volbaALFs==1 || volbaALFs==3
                                Wrpor=AWrpor*cosla+BWrpor*sinla;
                                clear AWrpor BWrpor
                                Wlambdapor=-AWlambdapor*sinla+BWlambdapor*cosla;
                                clear AWlambdapor BWlambdapor
                                Wfipor=AWfipor*cosla+BWfipor*sinla;
                                clear AWfipor BWfipor
                                Ur=AUr*cos(0*lambda');
                                clear AUr
                                Ufi=AUfi*cos(0*lambda');
                                clear AUfi
                            elseif volbaALFs==2
                                Wrpor=Wrpor*1e280;
                                Wlambdapor=Wlambdapor*1e280;
                                Wfipor=Wfipor*1e280;
                                Ur=Ur*1e280;
                                Ufi=Ufi*1e280;
                            end

                            Wrpor=bsxfun(@times,GM./r.^2,Wrpor);
                            Wrpor=bsxfun(@plus,GM./r.^2,Wrpor);
                            Wrpor=sqrt(bsxfun(@plus,-Wrpor,omegaEl^2.*r.*(cos(fiG).^2)).^2);                    
                            Wrpor=Wrpor(:);

                            Wlambdapor=bsxfun(@times,GM./r,Wlambdapor);
                            Wlambdapor=bsxfun(@times,Wlambdapor,1./(r.*cos(fiG)));
                            Wlambdapor=Wlambdapor(:);

                            Wfipor=bsxfun(@times,GM./r,Wfipor);
                            Wfipor=bsxfun(@plus,Wfipor,-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                            Wfipor=bsxfun(@times,Wfipor,1./r);
                            Wfipor=Wfipor(:);

                            Ur=bsxfun(@times,GM./r.^2,Ur);
                            Ur=bsxfun(@plus,GM./r.^2,Ur);
                            Ur=sqrt(bsxfun(@plus,-Ur,omegaEl^2.*r.*(cos(fiG).^2)).^2);                    
                            Ur=Ur(:);

                            Ufi=bsxfun(@times,GM./r,Ufi);
                            Ufi=bsxfun(@plus,Ufi,-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                            Ufi=bsxfun(@times,Ufi,1./r);
                            Ufi=Ufi(:);

                            Pg=(sqrt(Wrpor.^2+Wlambdapor.^2+Wfipor.^2)-sqrt(Ur.^2+Ufi.^2))*10^5;

                            clear Wrpor Wlambdapor Wfipor Ur Ufi
                        elseif volbapar(i)==21 %Gravity disturbance sa
                            if volbaALFs==1 || volbaALFs==3
                                porucha_sa=Aporucha_sa*cosla+Bporucha_sa*sinla;
                                clear Aporucha_sa Bporucha_sa
                            elseif volbaALFs==2
                                porucha_sa=porucha_sa*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^2,porucha_sa)*10^5;
                            Pg=Pg(:);

                            clear porucha_sa 
                        elseif volbapar(i)==22 %Height anomaly ell
                            if volbaALFs==1 || volbaALFs==3
                                zetaEl=AzetaEl*cosla+BzetaEl*sinla;
                                clear AzetaEl BzetaEl
                            elseif volbaALFs==2
                                zetaEl=zetaEl*1e280;
                            end

                            Pg=bsxfun(@times,GM./(r.*gamaP),zetaEl);
                            Pg=Pg(:);

                            clear zetaEl 
                        elseif volbapar(i)==23 %Height anomaly
                            
                            clear HC HCm HS HSm LmH Lmdg
                            
                            if volbaALFs==1 || volbaALFs==3
                                zeta_H=AH_zeta*cosla+BH_zeta*sinla;
                                clear AH_zeta BH_zeta
                                zeta_N1c=AN1c_zeta*cosla+BN1c_zeta*sinla;
                                clear AN1c_zeta BN1c_zeta
                                zeta_dg=Azetadg*cosla+Bzetadg*sinla;
                                clear Azetadg Bzetadg
                                zeta_zetaEl=Azeta*cosla+Bzeta*sinla;
                                clear Azeta Bzeta  
                            elseif volbaALFs==2
                                zeta_H=zeta_H*1e280;
                                zeta_N1c=zeta_N1c*1e280;
                                zeta_dg=zeta_dg*1e280;
                                zeta_zetaEl=zeta_zetaEl*1e280;
                            end                                 
                            
                            zeta_N1c=bsxfun(@times,GM./(r.*gamaP),zeta_N1c);
                            zeta_H(zeta_H<0)=zeta_H(zeta_H<0)*0; %H is set to zero in the areas of oceans and seas
   
                            zeta_zetaEl=bsxfun(@times,GM./(r.*gamaP),zeta_zetaEl);
                            
                            zeta_dg=bsxfun(@times,GM./r.^2,zeta_dg);
                            
                            G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                            ro=2670; %Density of the crust
                            
                            zeta_N=(zeta_N1c-bsxfun(@times,(2*pi*G*ro*zeta_H.^2),1./gamaP));

                            Pg=zeta_zetaEl-bsxfun(@times,zeta_dg.*(zeta_H+zeta_N),1./gamaP);
                            Pg=Pg(:);
                            
                            clear zeta_N1c zeta_H zeta_zetaEl zeta_dg zeta_N
                        elseif volbapar(i)==24 %Second radial derivative of disturbing potential
                            if volbaALFs==1 || volbaALFs==3
                                T_rr=AT_rr*cosla+BT_rr*sinla;
                                clear AT_rr BT_rr
                            elseif volbaALFs==2
                                T_rr=T_rr*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^3,T_rr)*10^9;
                            Pg=Pg(:);

                            clear T_rr
                        elseif volbapar(i)==25 %Second radial derivative of gravity potential
                            if volbaALFs==1 || volbaALFs==3
                                Wrr=AWrr*cosla+BWrr*sinla;
                                clear AWrr BWrr
                            elseif volbaALFs==2
                                Wrr=Wrr*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^3,Wrr);
                            clear Wrr
                            if nmin0==1 %Zero degree term
                                Pg=bsxfun(@plus,2*GM./r.^3,Pg);
                            end
                            Pg=bsxfun(@plus,Pg,omegaEl^2.*cos(fiG).^2)*10^9;
                            Pg=Pg(:);
                        end

                        if i==1
                            P=Pg;
                            clear Pg
                        else
                            P=[P Pg];
                            if i==pocetpar
                                clear Pg
                            end
                        end
                    end

                    %Update of the progress bar
                    set(progressbar,'string','','fontsize',8); drawnow;
                    
                    clear sinla cosla r gamaP                                
                end                

                %% Computation of functionals in point-wise
                %==============================================================

                %Identification of point-wise approach                         
                if volbadiskcheck==1                           
                    %Ellipsoidal coordinates
                    fi=str2num(get(findobj('tag','fi'),'string'))';
                    lambda=str2num(get(findobj('tag','lambda'),'string'))';
                    h=str2num(get(findobj('tag','hdisk'),'string'))';
                end

                %Identification of load data approach
                if volbaloadcheck==1
                    %Import of data file containing ellipsoidal coordinates
                    loadname=get(findobj('tag','use'),'userdata');
                    loadadresar=get(findobj('tag','diskcheck'),'userdata');

                    if isempty(loadname) %Error message, if GGM file has not been imported
                        errordlg('Please input the data file containing coordinates of the computing points.',...
                            'Error in point type selection');
                        error('Please input the data file containing coordinates of the computing points.')
                    end
                    
                    Import=load([loadadresar,loadname],'-ascii');

                    %Ellipsoidal coordinates
                    fi=Import(:,1);
                    lambda=Import(:,2);
                    h=Import(:,3);

                    volbadiskcheck=1;
                end

                %=============================================================
                if volbadiskcheck==1                             

                    %Error message for inputted ellipsoidal coordinates
                    if isempty(fi) || isempty(lambda) || isempty(h)
                        errordlg('At least one ellipsoidal coordinate is empty or incorrectly entered.',...
                            'Error in point type selection');
                        error('At least one ellipsoidal coordinate is empty or incorrectly entered.');
                    end

                    if length(fi)~=length(lambda) || length(fi)~=length(h) || length(lambda)~=length(h)
                        errordlg('Ellipsoidal coordinates dimensions are not consistent.',...
                            'Error in point type selection')
                        error('Ellipsoidal coordinates dimensions are not consistent.')     
                    end               

                    if any(fi>90) || any(fi<-90) 
                        errordlg('Values of Latitude must be within the interval <-90�,90�>.',...
                            'Error in point type selection');
                        error('Values of Latitude must be within the interval <-90�,90�>.');
                    end
                    if any(lambda>360) || any(lambda<-180)
                        errordlg('Values of Longitude min must be within the interval <-180�,180�> or <0�,360�>.',...
                            'Error in point type selection');
                        error('Values of Longitude min must be within the interval <-180�,180�> or <0�,360�>.');
                    end

                    fi=deg2rad(fi(:));
                    lambda=deg2rad(lambda(:));               

                    if coord==1 %Entered spherical coordinates                       
                        %Spherical radius
                        r=h;

                        %Spherical latitude
                        fiG=fi;  
                        
                        %Transform spherical latitude into the ellipsoidal
                        %latitude
                        [X Y Z]=sph2cart(lambda,fiG,r);
                        [fi , ~, h]=ecef2geodetic(X,Y,Z,[aEl eEl]');

                        clear X Y Z
                    elseif coord==0 %Entered ellipsoidal coordinates
                        %Trasformation of (fi, lambda, h) into (X, Y, Z)
                        [X,Y,Z]=geodetic2ecef(fi,lambda,h,[aEl eEl]');  
                        r=sqrt(X.*X+Y.*Y+Z.*Z); %Radius

                        %Spherical latitude
                        fiG=atan(Z./sqrt(X.*X+Y.*Y)); 

                        clear X Y Z 
                    end

                    %Computation of the coefficients C2,0; C4,0; C6,0;
                    %C8,0; C10,0 of the selected ellipsoid
                    CEl=zeros(length(C),1);
                    for n=1:5
                        CEl(2*n==stupen & rad==0,1)=((-1)^n*(3*eEl^(2*n))/((2*n+1)*(2*n+3)*sqrt(4*n+1))*(1-n-5^(3/2)*n*CEl_20/eEl^2)).*(aEl./R).^(2*n).*(GMEl/GM);
                    end                                

                    if any(volbapar==11) || any(volbapar==12) || any(volbapar==13) || any(volbapar==14) || any(volbapar==15) || any(volbapar==16) || any(volbapar==17) || any(volbapar==18) || any(volbapar==20) || any(volbapar==25)
                        grav=1;
                    else
                        grav=0;
                    end

                    if  any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==5) || any(volbapar==6) || any(volbapar==7) || any(volbapar==8) || any(volbapar==9) || any(volbapar==10) || any(volbapar==19) || any(volbapar==21) || any(volbapar==22) || any(volbapar==23) || any(volbapar==24)
                        por=1;
                        deltaC=C-CEl;
                    else
                        por=0;
                    end

                    if any(volbapar==20)
                        normal=1;
                    else
                        normal=0;
                    end

                    clear GGM stupen rad
                    if normal==0
                        clear CEl
                    end                                           

                    %Initialization
                    eta=0; ksi=0; Teta=0; Tksi=0; T=0; T_rr=0; Trr=0; Trfi=0; Trl=0; Tfifi=0;
                    Tfil=0; Tll=0; N=0; V=0; Vrr=0; Vrfi=0; Vrl=0; Vfifi=0;
                    Vfil=0; Vll=0; g=0; g_sa=0; W=0; anomalia_sa=0; porucha=0; 
                    porucha_sa=0; zetaEl=0; zeta=0; Wrr=0; Wr=0; Wfi=0; Wlambda=0;
                    Ur=0; Ufi=0; N1c=0; N2c=0; H=0; zeta_N1c=0; zeta_H=0; zeta_dg=0;
                    etaH=0; ksiH=0; TetaH=0; TksiH=0; TH=0; T_rrH=0; TrrH=0; TrfiH=0; TrlH=0; TfifiH=0;
                    TfilH=0; TllH=0; NH=0; VH=0; VrrH=0; VrfiH=0; VrlH=0; VfifiH=0;
                    VfilH=0; VllH=0; gH=0; g_saH=0; WH=0; anomalia_saH=0; poruchaH=0; 
                    porucha_saH=0; zetaElH=0; zetaH=0; WrrH=0; WrH=0; WfiH=0; WlambdaH=0;
                    UrH=0; UfiH=0; N1cH=0; N2cH=0; HH=0; zeta_N1cH=0; zeta_HH=0; zeta_dgH=0;
                    Tzz=0; Txx=0; Tyy=0; Txy=0; Txz=0; Tyz=0;
                    TzzH=0; TxxH=0; TyyH=0; TxyH=0; TxzH=0; TyzH=0;
                    Vzz=0; Vxx=0; Vyy=0; Vxy=0; Vxz=0; Vyz=0;
                    VzzH=0; VxxH=0; VyyH=0; VxyH=0; VxzH=0; VyzH=0;
                    
                    %Indices of the spherical harmonic coefficients
                    if radenie==0
                        index=zeros(nmax-1,1);
                        index(1,1)=1;
                        for i=1:(nmax-2)
                            index(i+1,1)=index(i,1)+2+i;
                        end
                    elseif radenie==1
                        z=((nmaxGGM+1)*(nmaxGGM+2)-6)/2-(sum(nmaxGGM-(nmaxGGM:-1:(nmax-1)))-1);
                        k=z;
                    end

                    %Initialization of the matrices and vectors for the computation of fnALFs
                    Pnm=zeros(length(fi),nmax+1);
                    q=(R./r);
                    q2=(R./r).^2;
                    u=cos(fiG);
                    t=sin(fiG);
                    %Initialization for extended-range arithmetic approach
                    if volbaALFs==3
                        
                        bit=mexext; %Bit version of Matlab
                        bit=bit(end-1:end);
                        bit=str2double(bit);
                        
                        nmax23=nmax*2+3;
                        rr=zeros(nmax23,1); ri=rr;
                        dd=zeros(nmax,1); am=dd; bm=am;
                        
                        if bit==32
                            pm=am;
                        else
                            zz=zeros(length(fiG),1);
                        end
                        
                        ps1=zeros(length(fiG),nmax); ips1=ps1;

                        m1=1:nmax23;
                        rr(m1)=sqrt(m1);
                        ri(m1)=1./rr;
                        m2=1:nmax;
                        dd(m2)=rr(2*m2+3).*ri(2*m2+2);

                        IND=960;
                        BIG=2^IND;
                        BIGI=2^(-IND);
                        BIGS=2^(IND/2);
                        BIGSI=2^(-IND/2);
                        ROOT3=1.732050807568877;
                        x=ROOT3*u.*q;
                        ix=zeros(size(x));
                        ps1(:,1)=x;
                        ips1(:,1)=ix;
                        for m3=2:nmax
                            x=(dd(m3-1)*u).*x.*q;
                            y=abs(x);
                            iy=y>=BIGS;
                            if any(iy)
                                x(iy)=x(iy)*BIGI;
                                ix(iy)=ix(iy)+1;
                            end
                            iy=y<BIGSI;
                            if any(iy)
                                x(iy)=x(iy)*BIG;
                                ix(iy)=ix(iy)-1;
                            end
                            ps1(:,m3)=x;
                            ips1(:,m3)=ix;
                        end
                    end

                    %If geoid/quasigeoid is to be computed
                    if any(volbapar==10) || any(volbapar==23)
                        geoid=1;
                        if any(h~=0)
                            errordlg('In order to compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                                'Error in point type selection');
                            error('In order to compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                        end
                    else
                        geoid=0;
                    end

                    %Initialization of the matrices and vectors for the 
                    %computation of the first-order derivatives of fnALFs
                    if any(volbapar==3) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16)  || any(volbapar==20)
                        dALFs=1;
                        dPnm=zeros(length(fi),nmax+1);
                        qu=q./u;
                        tu=t./u;
                        
                        %Treatment of the dPnm singularity
                        singdPnm=fi==pi/2 | fi==-pi/2;
                    else
                        dALFs=0;
                    end   
                    
                    %Initialization of the matrices and vectors for the 
                    %computation of the second-order derivatives of fnALFs
                    if any(volbapar==6) || any(volbapar==12)
                        ddALFs=1;
                        ddPnm=zeros(length(fi),nmax+1);
                        
                        %Treatment of the ddPnm singularity
                        singddPnm=fi==pi/2 | fi==-pi/2;
                    else
                        ddALFs=0;
                    end   

                    %Status line
                    progressbar=findobj('tag','hlasky');

                    %% Summation over m
                    for m=nmax:-1:0

                        %Update of the progress bar
                        if rem(m,10)==0
                            set(progressbar,'string',...
                                sprintf('Progress: m = %5.0d',m),...
                                'fontsize',8); drawnow;
                        end

                        %Selection of the spherical harmonic coefficients of order m
                        %======================================================
                        if m<2 
                            mu=2;
                            if radenie==0                             
                                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                                   Cm=C(index+m);
                                end

                                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                                    deltaCm=deltaC(index+m);
                                end

                                if normal==1 %C's spherical harmonic coefficients for the functionals with the normal field
                                    if m==0                                   
                                        CElm=CEl(index+m);
                                    end
                                end

                                if geoid==1
                                    if m==1
                                        HCm=HC([3;index+m+3]); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                        HSm=HS([3;index+m+3]); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                    elseif m==0
                                        HCm=HC([1;2;index+m+3]); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                        HSm=HS([1;2;index+m+3]); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                    end
                                end
                                
                                Sm=S(index+m);
                            elseif radenie==1
                                z=z+1;

                                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                                    Cm=C(z:k);
                                end

                                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                                    deltaCm=deltaC(z:k);
                                end

                                if normal==1 %C's spherical harmonic coefficients for the functionals with the normal field
                                    if m==0
                                        CElm=CEl(z:k);
                                    end
                                end
                                
                                if geoid==1
                                    if m==1
                                        HCm=HC((z-1):k);
                                        HSm=HS((z-1):k);
                                    elseif m==0
                                        HCm=HC((z-1):k);
                                        HSm=HS((z-1):k);
                                    end
                                end
                                
                                Sm=S(z:k);

                                zT=z;
                                kT=k;
                                
                                z=z-nmaxGGM+m-1;
                                k=z+nmax-m; 
                            end
                        else
                            mu=m;
                            if radenie==0
                                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                                    Cm=C(index((m-1):end)+m);
                                end

                                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                                    deltaCm=deltaC(index((m-1):end)+m);
                                end

                                if geoid==1
                                    HCm=HC(index((m-1):end)+m+3); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                    HSm=HS(index((m-1):end)+m+3); % '+3', since DMR also contains the coefficiets 00, 10, 11
                                end
                                
                                Sm=S(index((m-1):end)+m);
                            elseif radenie==1
                                if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                                    Cm=C(z:k);
                                end

                                if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                                    deltaCm=deltaC(z:k);
                                end

                                if geoid==1
                                    HCm=HC(z:k);
                                    HSm=HS(z:k);
                                end
                                
                                Sm=S(z:k);
                                
                                zT=z;
                                kT=k;
                                
                                z=z-nmaxGGM+m-2;
                                k=z+nmax-m+1; 
                            end
                        end
                        %======================================================

                        %% Computation of modified fnALFs
                        if volbaALFs==1 %Standard forward column method
                            if m==0
                                Pnm(:,1)=1;
                            elseif m==1                   
                                Pnm(:,1)=sqrt(3)*u.*q;  
                            elseif m>1                            
                                i=2*(2:m);
                                i1=sqrt((i+ones(size(i)))./i);
                                Pnm(:,1)=u.^m*sqrt(3)*prod(i1).*q.^m;
                            end

                            if m==nmax
                            elseif m<=(nmax-1)
                                n=m+1;
                                anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                            end

                            if m<(nmax-1)
                                j=3;
                                for n=m+2:nmax                            
                                    anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                    bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                                    Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                                    j=j+1;
                                end
                            end

                        elseif volbaALFs==2 %Modified forward column method
                            if m==0
                                Pnm(:,1)=1e-280;
                            elseif m==1

                                Pnm(:,1)=sqrt(3)*q*1e-280;  
                            elseif m>1                            
                                i=2*(2:m);
                                i1=sqrt((i+ones(size(i)))./i);
                                Pnm(:,1)=sqrt(3)*prod(i1)*(q.^m)*1e-280;
                            end

                            if m==nmax
                            elseif m<=(nmax-1)
                                n=m+1;
                                anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                            end

                            if m<(nmax-1)
                                j=3;
                                for n=m+2:nmax                            
                                    anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                    bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                                    Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                                    j=j+1;
                                end
                            end

                        elseif volbaALFs==3 %Extended-range arithmetic 
                            if bit==32 %32 bit version of Matlab
                                am(m+1)=rr(2*m+3);
                                for n=m+2:nmax
                                    w=rr(2*n+1)*ri(n-m)*ri(n+m);
                                    am(n)=rr(2*n-1)*w;
                                    bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*w;
                                end

                                if m~=0
                                    for i=1:length(fiG)                                   
                                        x=ps1(i,m);
                                        ix=ips1(i,m);
                                        if(ix==0)
                                            pm(m)=x;
                                        elseif (ix<-1)    
                                            pm(m)=0;         
                                        elseif (ix<0)
                                            pm(m)=x*BIGI;
                                        else
                                            pm(m)=x*BIG;
                                        end

                                        if(m>=nmax)
                                            Pnm(i,1:(nmax-m+1))=pm(m:end);
                                            continue;
                                        end
                                        y=x;
                                        iy=ix;
                                        x=(am(m+1)*t(i)*q(i))*y;
                                        ix=iy;
                                        w=abs(x);

                                        if(w>=BIGS)
                                            x=x*BIGI;
                                            ix=ix+1;
                                        elseif (w<BIGSI)
                                            x=x*BIG;
                                            ix=ix-1;
                                        end
                                        if(ix==0)
                                            pm(m+1)=x;
                                        elseif (ix<-1)    
                                            pm(m+1)=0.;        
                                        elseif (ix<0)
                                            pm(m+1)=x*BIGI;
                                        else
                                            pm(m+1)=x*BIG;
                                        end

                                        for n=m+2:nmax                                       
                                            id=ix-iy;
                                            if(id==0)
                                                zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*y;
                                                iz=ix;
                                            elseif (id==1)
                                                zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*(y*BIGI);
                                                iz=ix;
                                            elseif (id==-1)
                                                zz=(am(n)*t(i)*q(i))*(x*BIGI)-bm(n)*q2(i)*y;
                                                iz=iy;
                                            elseif (id>1)
                                                zz=(am(n)*t(i)*q(i))*x;
                                                iz=ix;
                                            else
                                                zz=-bm(n)*q2(i)*y;
                                                iz=iy;
                                            end

                                            w=abs(zz);
                                            if(w>=BIGS)
                                                zz=zz*BIGI;
                                                iz=iz+1;
                                            elseif (w<BIGSI)
                                                zz=zz*BIG;
                                                iz=iz-1;
                                            end

                                            if(iz==0)
                                                pm(n)=zz;
                                            elseif (iz<-1)   
                                                pm(n)=0.;          
                                            elseif (iz<0)
                                                pm(n)=zz*BIGI;
                                            else
                                                pm(n)=zz*BIG;
                                            end

                                            y=x;
                                            iy=ix;
                                            x=zz;
                                            ix=iz; 
                                        end

                                        Pnm(i,1:(nmax-m+1))=pm(m:end);         
                                    end

                                elseif m==0
                                    Pnm(:,1)=1;
                                    Pnm(:,2)=sqrt(3)*t.*q;

                                    for i=2:nmax
                                        Pnm(:,i+1)=Pnm(:,i).*sqrt((2*i+1)*(2*i-1))./i.*t.*q-q2.*Pnm(:,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                                    end

                                    clear rr ri dd am bm pm ps1 ips1 m1 m2 ...
                                        dd ix x y iy w iz z
                                end
                            else %64 bit version of Matlab
                                am(m+1)=rr(2*m+3);
                                for n=m+2:nmax
                                    w=rr(2*n+1)*ri(n-m)*ri(n+m);
                                    am(n)=rr(2*n-1)*w;
                                    bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*w;
                                end
                                                              
                                if m==0 %Zonal modified fnALFs
                                    Pnm(:,1)=1;
                                    Pnm(:,2)=sqrt(3)*t.*q;

                                    for i=2:nmax
                                        Pnm(:,i+1)=Pnm(:,i).*sqrt((2*i+1)*(2*i-1))./i.*t.*q-q2.*Pnm(:,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                                    end

                                    clear rr ri am bm ps1 ips1 m1 m2 dd ...
                                        ix x y iy w iz zz pmx pm0 pmxBIGI ...
                                        pmxBIG w wBIGS wBIGSI pm1x pm10 ...
                                        pm1xBIGI pm1xBIG id id0 id1 id_1 ...
                                        idv1 idm1 iz0 izm_1 izm0 izv0

                                elseif m~=0 %Non-zonal modified fnALFs
                                   x=ps1(:,m);
                                   ix=ips1(:,m);

                                   pmx=ix==0;
                                   pm0=ix<-1;
                                   pmxBIGI=(ix>=-1 & ix<0);
                                   pmxBIG=ix>0;

                                   Pnm(pmx,1)=x(pmx);
                                   Pnm(pm0,1)=0;
                                   Pnm(pmxBIGI,1)=x(pmxBIGI)*BIGI;
                                   Pnm(pmxBIG,1)=x(pmxBIG)*BIG; 
                                   
                                   if m<nmax
                                       y=x;
                                       iy=ix;

                                       x=(am(m+1).*t).*y.*q;
                                       ix=iy;
                                       w=abs(x);

                                       wBIGS=w>=BIGS;
                                       x(wBIGS)=x(wBIGS)*BIGI;
                                       ix(wBIGS)=ix(wBIGS)+1;

                                       wBIGSI=w<BIGSI;
                                       x(wBIGSI)=x(wBIGSI)*BIG;
                                       ix(wBIGSI)=ix(wBIGSI)-1;   

                                       pm1x=ix==0;
                                       pm10=ix<-1;
                                       pm1xBIGI=(ix>=-1 & ix<0);
                                       pm1xBIG=ix>0;

                                       Pnm(pm1x,2)=x(pm1x);
                                       Pnm(pm10,2)=0;
                                       Pnm(pm1xBIGI,2)=x(pm1xBIGI)*BIGI;
                                       Pnm(pm1xBIG,2)=x(pm1xBIG)*BIG;

                                       for n=m+2:nmax
                                           id=ix-iy;

                                           id0=id==0;
                                           id1=id==1;
                                           id_1=id==-1;
                                           idv1=id>1;
                                           idm1=id<-1;

                                           zz(id0)=(am(n).*t(id0).*q(id0)).*x(id0)-bm(n).*y(id0).*q2(id0);
                                           iz(id0,1)=ix(id0);

                                           zz(id1)=(am(n).*t(id1).*q(id1)).*x(id1)-bm(n).*(y(id1).*q2(id1).*BIGI);
                                           iz(id1)=ix(id1);

                                           zz(id_1)=(am(n).*t(id_1).*q(id_1)).*(x(id_1).*BIGI)-bm(n).*q2(id_1).*y(id_1);
                                           iz(id_1)=iy(id_1);

                                           zz(idv1)=(am(n).*t(idv1).*q(idv1)).*x(idv1);
                                           iz(idv1)=ix(idv1);

                                           zz(idm1)=-bm(n).*y(idm1).*q2(idm1);
                                           iz(idm1)=iy(idm1);

                                           w=abs(zz);

                                           wBIGS=w>=BIGS;
                                           zz(wBIGS)=zz(wBIGS)*BIGI;
                                           iz(wBIGS)=iz(wBIGS)+1;

                                           wBIGSI=w<BIGSI;
                                           zz(wBIGSI)=zz(wBIGSI)*BIG;
                                           iz(wBIGSI)=iz(wBIGSI)-1;

                                           iz0=iz==0;
                                           izm_1=iz<-1;
                                           izm0=(iz>=-1 & iz<0);
                                           izv0=iz>0;

                                           Pnm(iz0,n-m+1)=zz(iz0);
                                           Pnm(izm_1,n-m+1)=0;
                                           Pnm(izm0,n-m+1)=zz(izm0)*BIGI;
                                           Pnm(izv0,n-m+1)=zz(izv0)*BIG;                           
                                                         
                                           y=x;
                                           iy=ix;
                                           x=zz;
                                           ix=iz;
                                       end   
                                    end
                                end
                            end
                        end

                        %If nmin~=2
                        %======================================================
                        if nmin~=2
                            if LNOFnmin==0
                                if m<2
                                    if por==1
                                        deltaCm(1:(nmin-2))=0;
                                    end

                                    if grav==1
                                        Cm(1:(nmin-2))=0;
                                    end

                                    if geoid==1
                                        HCm(1:(nmin-m))=0;
                                        HSm(1:(nmin-m))=0;
                                    end

                                    if normal==1
                                        if m==0
                                            CElm(1:(nmin-2))=0;
                                        end
                                    end

                                    Sm(1:(nmin-2))=0;
                                else
                                    if por==1
                                        deltaCm(1:(nmin-m))=0;
                                    end

                                    if grav==1
                                        Cm(1:(nmin-m))=0;
                                    end

                                    if geoid==1
                                        HCm(1:(nmin-m))=0;
                                        HSm(1:(nmin-m))=0;
                                    end

                                    Sm(1:(nmin-m))=0;
                                end
                            elseif LNOFnmin==1
                                Pnm(:,1:(nmin-m))=0;
                            end
                        end
                        %======================================================
                        
                        %% Computation of the first-order derivatives of modified fnALFs
                        if dALFs==1  
                            if volbaALFs==1 || volbaALFs==3
                                enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                                if m==0 %Zonal modified dALFs
                                    dPnm(:,1)=0.*u;
                                    dPnm(:,2)=sqrt(3)*u.*q;
                                    dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                                elseif m==nmax %Sectorial modified dALFs
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1);
                                else
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                    dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne ALFs
                                end

                            elseif volbaALFs==2
                                enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                                if m==0 %Zonal modified dALFs
                                    dPnm(:,1)=0.*u*1e-280;
                                    dPnm(:,2)=sqrt(3)*u.*q*1e-280;
                                    dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                                elseif m==nmax
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Sectorial modified dALFs
                                else
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                    dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne dALFs
                                end
                            end                   

                            %Treatment of the dALFs singularity
                            dPnm(singdPnm,:)=0;
                            
                            if ddALFs==1 %If the second-order derivatives of the modified fnALFs are to be computed
                                if m==0 %Zonal modified ddALFs
                                    ddPnm=bsxfun(@times,tu,dPnm)-bsxfun(@times,(0:nmax).*((0:nmax)+1),Pnm);
                                else
                                    
                                    ddPnm(:,1:end-m)=bsxfun(@times,tu,dPnm(:,1:end-m))+bsxfun(@times,m^2./u.^2,Pnm(:,1:end-m))-bsxfun(@times,(m:nmax).*((m:nmax)+1),Pnm(:,1:end-m));
                                end
                                                                
                                %Treatment of the ddALFs singularity
                                ddPnm(singddPnm,:)=0;
                            end
                            
                        end                                     

                        cosla=cos(m*lambda);
                        sinla=sin(m*lambda);

                        if any(volbapar==8) || any(volbapar==9) || any(volbapar==14) || any(volbapar==15)
                            if m<2
                                coslaplus2=cos((m+2)*lambda);
                                sinlaplus2=sin((m+2)*lambda);
                                
                                coslaminus2=cos(0*lambda);
                                sinlaminus2=sin(0*lambda);
                            elseif m>nmax-2
                                coslaplus2=cos(0*lambda);
                                sinlaplus2=sin(0*lambda);
                                
                                coslaminus2=cos((m-2)*lambda);
                                sinlaminus2=sin((m-2)*lambda);
                            else
                                coslaplus2=cos((m+2)*lambda);
                                sinlaplus2=sin((m+2)*lambda);
                                
                                coslaminus2=cos((m-2)*lambda);
                                sinlaminus2=sin((m-2)*lambda);
                            end
                        end
                        
                        if any(volbapar==9) || any(volbapar==15)
                            if m<1
                                coslaplus1=cos((m+1)*lambda);
                                sinlaplus1=sin((m+1)*lambda);
                                
                                coslaminus1=cos(0*lambda);
                                sinlaminus1=sin(0*lambda);
                            elseif m>nmax-1
                                coslaplus1=cos(0*lambda);
                                sinlaplus1=sin(0*lambda);
                                
                                coslaminus1=cos((m-1)*lambda);
                                sinlaminus1=sin((m-1)*lambda);
                            else
                                coslaplus1=cos((m+1)*lambda);
                                sinlaplus1=sin((m+1)*lambda);
                                
                                coslaminus1=cos((m-1)*lambda);
                                sinlaminus1=sin((m-1)*lambda);
                            end
                        end
                        
                        %% Loop for 1:NF (number of computing functionals)                        
                        for i=1:pocetpar 

                            %Summation over n
                            if volbapar(i)==1         
                            elseif volbapar(i)==2 %Deflection of the vertical eta                                
                                 
                                if m<2  
                                    for n=mu:nmax
                                        eta=eta+(-deltaCm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*(m*Pnm(:,n-m+1));                                    
                                    end
                                else
                                    for n=mu:nmax
                                        eta=eta+(-deltaCm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*(m*Pnm(:,n-mu+1));
                                    end
                                end

                                if volbaALFs==2
                                    etaH=etaH.*u+eta;
                                    eta=0;
                                end

                            elseif volbapar(i)==3 %Deflection of the vertical xi
                             
                                if m<2  
                                    for n=mu:nmax
                                        ksi=ksi+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*dPnm(:,n-m+1);                                    
                                    end
                                else
                                    for n=mu:nmax
                                        ksi=ksi+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*dPnm(:,n-mu+1);
                                    end
                                end
                                
                                if volbaALFs==2
                                    ksiH=ksiH.*u+ksi;
                                    ksi=0;
                                end

                            elseif volbapar(i)==4 %Deflection of the vertical Theta
                              
                                if m<2     
                                    for n=mu:nmax
                                        Teta=Teta+(-deltaCm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*(m*Pnm(:,n-m+1));
                                        Tksi=Tksi+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*dPnm(:,n-m+1); 
                                    end
                                else
                                    for n=mu:nmax
                                        Teta=Teta+(-deltaCm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*(m*Pnm(:,n-mu+1));
                                        Tksi=Tksi+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*dPnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    TetaH=TetaH.*u+Teta;
                                    Teta=0;
                                    TksiH=TksiH.*u+Tksi;
                                    Tksi=0;
                                end

                            elseif volbapar(i)==5 %Disturbing potential
                             
                                if m<2 
                                    for n=mu:nmax
                                        T=T+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);  
                                    end
                                else
                                    for n=mu:nmax
                                        T=T+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    TH=TH.*u+T;
                                    T=0;
                                end

                            elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                             
                                if m<2  
                                    for n=mu:nmax
                                        deltaCmcosla=deltaCm(n-mu+1).*cosla;
                                        Smsinla=Sm(n-mu+1).*sinla;

                                        Trr=Trr+(n+1).*(n+2).*(deltaCmcosla+Smsinla).*Pnm(:,n-m+1);                                    
                                        Tfifi=Tfifi+(deltaCmcosla+Smsinla).*ddPnm(:,n-m+1);
                                        Tll=Tll+(deltaCmcosla+Smsinla).*m.^2.*Pnm(:,n-m+1);
                                    end
                                else
                                    for n=mu:nmax
                                        deltaCmcosla=deltaCm(n-mu+1).*cosla;
                                        Smsinla=Sm(n-mu+1).*sinla;

                                        Trr=Trr+(n+1).*(n+2).*(deltaCmcosla+Smsinla).*Pnm(:,n-mu+1);
                                        Tfifi=Tfifi+(deltaCmcosla+Smsinla).*ddPnm(:,n-mu+1);
                                        Tll=Tll+(deltaCmcosla+Smsinla).*m.^2.*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    TrrH=TrrH.*u+Trr;
                                    Trr=0;
                                    TfifiH=TfifiH.*u+Tfifi;
                                    Tfifi=0;
                                    TllH=TllH.*u+Tll;
                                    Tll=0;
                                end
                                
                            elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                                                              
                                if m<2   
                                    for n=mu:nmax
                                        Smcosla=Sm(n-mu+1).*cosla;
                                        deltaCmsinla=deltaCm(n-mu+1).*sinla;

                                        Trfi=Trfi+(n+1).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*dPnm(:,n-m+1);                                    
                                        Trl=Trl+(n+1).*(Smcosla-deltaCmsinla).*(m*Pnm(:,n-m+1));
                                        Tfil=Tfil+(Smcosla-deltaCmsinla).*(m*dPnm(:,n-m+1));
                                    end
                                else
                                    for n=mu:nmax
                                        Smcosla=Sm(n-mu+1).*cosla;
                                        deltaCmsinla=deltaCm(n-mu+1).*sinla;

                                        Trfi=Trfi+(n+1).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*dPnm(:,n-mu+1);
                                        Trl=Trl+(n+1).*(Smcosla-deltaCmsinla).*(m*Pnm(:,n-mu+1));
                                        Tfil=Tfil+(Smcosla-deltaCmsinla).*(m*dPnm(:,n-mu+1));
                                    end
                                end

                                if volbaALFs==2
                                    TrfiH=TrfiH.*u+Trfi;
                                    Trfi=0;
                                    TrlH=TrlH.*u+Trl;
                                    Trl=0;
                                    TfilH=TfilH.*u+Tfil;
                                    Tfil=0;
                                end
                                
                            elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                                
                                %Coefficients of the ALFs
                                if m==0
                                    anm=sqrt(2)/4.*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)-(m+2)+2);
                                    bnm=((2:nmax)+m+1).*((2:nmax)+m+2)./2./(m+1);
                                    cnm=zeros(nmax-1,1);
                                elseif m==1
                                    anm=1/4.*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)-(m+2)+2);
                                    bnm=((2:nmax)+m+1).*((2:nmax)+m+2)./2./(m+1);
                                    cnm=zeros(nmax-1,1);
                                elseif m==2
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    cnm=sqrt(2)/4.*sqrt((2:nmax).^2-(m-2+1).^2).*sqrt((2:nmax)-(m-2)).*sqrt((2:nmax)+m-2+2);
                                else
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    cnm=1/4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                end
                                  
                                %Spherical harmonic coefficients of orders m+2 and m-2
                                if radenie==0
                                    if m<2
                                        deltaCplus2=deltaC(index+m+2);
                                        Splus2=S(index+m+2); 
                                        
                                        deltaCminus2=zeros(nmax-1,1);
                                        Sminus2=zeros(nmax-1,1); 
                                    elseif m>nmax-2
                                        deltaCplus2=zeros(nmax-m+1,1);
                                        Splus2=zeros(nmax-m+1,1);
                                        
                                        deltaCminus2=deltaC(index((m-1):end)+m-2);
                                        Sminus2=S(index((m-1):end)+m-2); 
                                    else
                                        deltaCplus2=deltaC(index((m-1):end)+m+2);
                                        Splus2=S(index((m-1):end)+m+2);
                                        
                                        deltaCminus2=deltaC(index((m-1):end)+m-2);
                                        Sminus2=S(index((m-1):end)+m-2); 
                                    end
                                elseif radenie==1
                                    if m<2
                                        deltaCplus2=[zeros(m,1);deltaC((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))];
                                        Splus2=[zeros(m,1);S((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))]; 
                                        
                                        deltaCminus2=zeros(nmax-1,1);
                                        Sminus2=zeros(nmax-1,1); 
                                    elseif m>nmax-2
                                        if m==2
                                            deltaCminus2=deltaC(1:(nmax-1));
                                            Sminus2=deltaC(1:(nmax-1))*0;
                                        else
                                            deltaCminus2=deltaC((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            Sminus2=S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                        
                                        deltaCplus2=zeros(nmax-m+1,1);
                                        Splus2=zeros(nmax-m+1,1);
                                    else
                                        if m==2
                                            deltaCminus2=deltaC(1:(nmax-1));
                                            Sminus2=deltaC(1:(nmax-1))*0;
                                        else
                                            deltaCminus2=deltaC((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            Sminus2=S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                        
                                        deltaCplus2=[0;0;deltaC((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                        Splus2=[0;0;S((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                    end
                                end
                                
                                %Summation over m                    
                                if m<2   
                                    for n=mu:nmax
                                        %Tzz
                                        Tzz=Tzz+(n+1).*(n+2).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);    

                                        %Txx    
                                        Txx1=(deltaCplus2(n-mu+1).*coslaplus2+Splus2(n-mu+1).*sinlaplus2).*anm(n-mu+1);
                                        Txx2_vnutro=deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla;
                                        Txx2=Txx2_vnutro.*(bnm(n-mu+1)-(n+1)*(n+2));
                                        Txx3=(deltaCminus2(n-mu+1).*coslaminus2+Sminus2(n-mu+1).*sinlaminus2).*cnm(n-mu+1);
                                        Txx=Txx+(Txx1+Txx2+Txx3).*Pnm(:,n-m+1);  

                                        %Tyy
                                        Tyy2=Txx2_vnutro.*bnm(n-mu+1);
                                        Tyy=Tyy+(Txx1+Tyy2+Txx3).*Pnm(:,n-m+1); 
                                    end
                                else
                                    for n=mu:nmax
                                        %Tzz
                                        Tzz=Tzz+(n+1).*(n+2).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);

                                        %Txx
                                        Txx1=(deltaCplus2(n-mu+1).*coslaplus2+Splus2(n-mu+1).*sinlaplus2).*anm(n-mu+1);
                                        Txx2_vnutro=deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla;
                                        Txx2=Txx2_vnutro.*(bnm(n-mu+1)-(n+1)*(n+2));
                                        Txx3=(deltaCminus2(n-mu+1).*coslaminus2+Sminus2(n-mu+1).*sinlaminus2).*cnm(n-mu+1);
                                        Txx=Txx+(Txx1+Txx2+Txx3).*Pnm(:,n-mu+1);

                                        %Tyy
                                        Tyy2=Txx2_vnutro.*bnm(n-mu+1);
                                        Tyy=Tyy+(Txx1+Tyy2+Txx3).*Pnm(:,n-mu+1); 
                                    end
                                end

                                if volbaALFs==2
                                    TzzH=TzzH.*u+Tzz;
                                    Tzz=0;
                                    TxxH=TxxH.*u+Txx;
                                    Txx=0;
                                    TyyH=TyyH.*u+Tyy;
                                    Tyy=0;
                                end
                                
                            elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                                                               
                                %Coefficients of the ALFs
                                if m==0
                                    dnm=1/4.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt(2).*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)+m+2-2);
                                    gnm=zeros(nmax-1,1);
                                    hnm=gnm;
                                    
                                    betanm=((2:nmax)+2)./2.*sqrt(1+ones(1,nmax-1)).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)-(m+1)+1);
                                    gamanm=gnm;
                                    
                                    minm=((2:nmax)+2)./2.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt(2).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)+m+1-1);
                                    ninm=gnm;
                                elseif m==1
                                    dnm=1/4.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)+m+2-2);
                                    gnm=-1/4*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax)+1).*sqrt((2:nmax)-1).*((2:nmax)+2);
                                    hnm=zeros(nmax-1,1);
                                    
                                    betanm=((2:nmax)+2)./2.*sqrt((2:nmax)+m+1).*sqrt((2:nmax)-(m+1)+1);
                                    gamanm=-((2:nmax)+2).*sqrt((2:nmax).*((2:nmax)+1)./2);
                                    
                                    minm=((2:nmax)+2)./2.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)+m+1-1);
                                    ninm=gnm;
                                elseif m==2
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=zeros(nmax-1,1);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                elseif m==3
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-3).*sqrt((m:nmax)-2).*sqrt((m:nmax)-1).*sqrt((m:nmax)+2);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                else
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)-(m-2)-2);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                end
                                
                                %Spherical harmonic coefficients of orders m+2, m-2, m+1 and m-1
                                if radenie==0
                                    if m<2
                                        deltaCplus2=deltaC(index+m+2);
                                        Splus2=S(index+m+2); 
                                        
                                        deltaCminus2=zeros(nmax-1,1);
                                        Sminus2=zeros(nmax-1,1); 
                                        
                                        deltaCplus1=deltaC(index+m+1);
                                        Splus1=S(index+m+1);
                                        
                                        if m==0
                                            deltaCminus1=zeros(nmax-1,1);
                                            Sminus1=zeros(nmax-1,1);
                                        else
                                            deltaCminus1=deltaC(index+m-1);
                                            Sminus1=S(index+m-1);
                                        end
                                    elseif m>nmax-2
                                        deltaCplus2=zeros(nmax-m+1,1);
                                        Splus2=zeros(nmax-m+1,1);
                                        
                                        deltaCminus2=deltaC(index((m-1):end)+m-2);
                                        Sminus2=S(index((m-1):end)+m-2); 
                                        
                                        if m==nmax
                                            deltaCplus1=zeros(nmax-m+1,1);
                                            Splus1=zeros(nmax-m+1,1);
                                        else
                                            deltaCplus1=deltaC(index((m-1):end)+m+1);
                                            Splus1=S(index((m-1):end)+m+1);
                                        end
                                        
                                        deltaCminus1=deltaC(index((m-1):end)+m-1);
                                        Sminus1=S(index((m-1):end)+m-1); 
                                    else
                                        deltaCplus2=deltaC(index((m-1):end)+m+2);
                                        Splus2=S(index((m-1):end)+m+2);
                                        
                                        deltaCminus2=deltaC(index((m-1):end)+m-2);
                                        Sminus2=S(index((m-1):end)+m-2); 
                                        
                                        deltaCplus1=deltaC(index((m-1):end)+m+1);
                                        Splus1=S(index((m-1):end)+m+1);
                                        
                                        deltaCminus1=deltaC(index((m-1):end)+m-1);
                                        Sminus1=S(index((m-1):end)+m-1); 
                                    end
                                elseif radenie==1
                                    if m<2
                                        deltaCplus2=[zeros(m,1);deltaC((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))];
                                        Splus2=[zeros(m,1);S((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))]; 
                                        
                                        deltaCminus2=zeros(nmax-1,1);
                                        Sminus2=zeros(nmax-1,1); 
                                        
                                        deltaCplus1=deltaC((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)));
                                        Splus1=S((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)));
                                        
                                        if m==1
                                            deltaCminus1=deltaC(1:(nmax-1));
                                            Sminus1=deltaC(1:(nmax-1))*0;
                                        else                                            
                                            deltaCminus1=zeros(nmax-1,1);
                                            Sminus1=zeros(nmax-1,1);
                                        end
                                    elseif m>nmax-2
                                        if m==2
                                            deltaCminus2=deltaC(1:(nmax-1));
                                            Sminus2=deltaC(1:(nmax-1))*0;
                                        else
                                            deltaCminus2=deltaC((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            Sminus2=S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                        
                                        deltaCplus2=zeros(nmax-m+1,1);
                                        Splus2=zeros(nmax-m+1,1);
                                        
                                        if m==nmax
                                            deltaCplus1=zeros(nmax-m+1,1);
                                            Splus1=zeros(nmax-m+1,1);
                                        else
                                            deltaCplus1=[0;deltaC((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                            Splus1=[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                        end
                                        
                                        deltaCminus1=deltaC((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                        Sminus1=S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))); 
                                    else
                                        if m==2
                                            deltaCminus2=deltaC(1:(nmax-1));
                                            Sminus2=deltaC(1:(nmax-1))*0;
                                        else
                                            deltaCminus2=deltaC((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            Sminus2=S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                        
                                        deltaCplus2=[0;0;deltaC((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                        Splus2=[0;0;S((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                        
                                        deltaCplus1=[0;deltaC((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                        Splus1=[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                        
                                        deltaCminus1=deltaC((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                        Sminus1=S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                    end
                                end
                                
                                PTnm=[zeros(length(fi),1) Pnm];
                                
                                %Summation over n                          
                                if m<2  
                                    for n=mu:nmax
                                        %Txy    
                                        Txy1=(-deltaCplus2(n-mu+1).*sinlaplus2+Splus2(n-mu+1).*coslaplus2).*dnm(n-mu+1);
                                        Txy2=(-deltaCm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*gnm(n-mu+1);
                                        Txy3=(-deltaCminus2(n-mu+1).*sinlaminus2+Sminus2(n-mu+1).*coslaminus2).*hnm(n-mu+1);
                                        Txy=Txy+(Txy1+Txy2+Txy3).*PTnm(:,n-m+1).*q;  

                                        %Txz
                                        Txz1=(deltaCplus1(n-mu+1).*coslaplus1+Splus1(n-mu+1).*sinlaplus1).*betanm(n-mu+1);
                                        Txz2=(deltaCminus1(n-mu+1).*coslaminus1+Sminus1(n-mu+1).*sinlaminus1).*gamanm(n-mu+1);
                                        Txz=Txz+(Txz1+Txz2).*Pnm(:,n-m+1);

                                        %Tyz
                                        Tyz1=(-deltaCplus1(n-mu+1).*sinlaplus1+Splus1(n-mu+1).*coslaplus1).*minm(n-mu+1);
                                        Tyz2=(-deltaCminus1(n-mu+1).*sinlaminus1+Sminus1(n-mu+1).*coslaminus1).*ninm(n-mu+1);
                                        Tyz=Tyz+(Tyz1+Tyz2).*PTnm(:,n-m+1).*q;  
                                    end
                                else
                                    for n=mu:nmax
                                        %Txy
                                        Txy1=(-deltaCplus2(n-mu+1).*sinlaplus2+Splus2(n-mu+1).*coslaplus2).*dnm(n-mu+1);
                                        Txy2=(-deltaCm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*gnm(n-mu+1);
                                        Txy3=(-deltaCminus2(n-mu+1).*sinlaminus2+Sminus2(n-mu+1).*coslaminus2).*hnm(n-mu+1);
                                        Txy=Txy+(Txy1+Txy2+Txy3).*PTnm(:,n-mu+1).*q;

                                        %Txz
                                        Txz1=(deltaCplus1(n-mu+1).*coslaplus1+Splus1(n-mu+1).*sinlaplus1).*betanm(n-mu+1);
                                        Txz2=(deltaCminus1(n-mu+1).*coslaminus1+Sminus1(n-mu+1).*sinlaminus1).*gamanm(n-mu+1);
                                        Txz=Txz+(Txz1+Txz2).*Pnm(:,n-mu+1);

                                        %Tyz
                                        Tyz1=(-deltaCplus1(n-mu+1).*sinlaplus1+Splus1(n-mu+1).*coslaplus1).*minm(n-mu+1);
                                        Tyz2=(-deltaCminus1(n-mu+1).*sinlaminus1+Sminus1(n-mu+1).*coslaminus1).*ninm(n-mu+1);
                                        Tyz=Tyz+(Tyz1+Tyz2).*PTnm(:,n-mu+1).*q;  
                                    end
                                end

                                if volbaALFs==2
                                    TxyH=TxyH.*u+Txy;
                                    Txy=0;
                                    TxzH=TxzH.*u+Txz;
                                    Txz=0;
                                    TyzH=TyzH.*u+Tyz;
                                    Tyz=0;
                                end
                                
                            elseif volbapar(i)==10 %Geoid undulation
        
                                if m<2 %Spherical harmonics of geoid undulation of orders 0 a 1 
                                    
                                    % When computing H, there is no
                                    % dumping factor (R./r).^n,
                                    % therefore the matrix Pnm has to be
                                    % devided by 1./((R./r).^n), since
                                    % Pnm is the matrix of the MODIFFIED
                                    % fnALFS
                                    for nn=m:1
                                        H=H+1./(R./r).^nn.*(HCm(nn-m+1).*cosla+HSm(nn-m+1).*sinla).*Pnm(:,nn-m+1); 
                                    end
                                end
                                
                                                                
                                if m<2 
                                    for n=mu:nmax
                                        N1c=N1c+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);  
                                        
                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFFIED
                                        % fnALFS
                                        H=H+1./(R./r).^n.*(HCm(n-mu+3-m).*cosla+HSm(n-mu+3-m).*sinla).*Pnm(:,n-m+1); 
                                    end
                                else
                                    for n=mu:nmax
                                        N1c=N1c+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                        
                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFFIED
                                        % fnALFS
                                        H=H+1./(R./r).^n.*(HCm(n-mu+1).*cosla+HSm(n-mu+1).*sinla).*Pnm(:,n-mu+1); 
                                    end
                                end
                                    
                                if volbaALFs==2
                                    N1cH=N1cH.*u+N1c;
                                    N1c=0;
                                    HH=HH.*u+H;
                                    H=0;
                                end
                             
                                                        
                            elseif volbapar(i)==11 %Gravitational potential
                               
                                if m<2       
                                    for n=mu:nmax 
                                        V=V+(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1); 
                                    end
                                else
                                    for n=mu:nmax 
                                        V=V+(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                    end
                                end
                                
                                if volbaALFs==2    
                                    VH=VH.*u+V;
                                    V=0;
                                end

                            elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                                
                                if m<2  
                                    for n=mu:nmax
                                        Cmcosla=Cm(n-mu+1).*cosla;
                                        Smsinla=Sm(n-mu+1).*sinla;

                                        Vrr=Vrr+(n+1).*(n+2).*(Cmcosla+Smsinla).*Pnm(:,n-m+1);                                    
                                        Vfifi=Vfifi+(Cmcosla+Smsinla).*ddPnm(:,n-m+1);
                                        Vll=Vll+(Cmcosla+Smsinla).*m.^2.*Pnm(:,n-m+1);
                                    end
                                else
                                    for n=mu:nmax
                                        Cmcosla=Cm(n-mu+1).*cosla;
                                        Smsinla=Sm(n-mu+1).*sinla;

                                        Vrr=Vrr+(n+1).*(n+2).*(Cmcosla+Smsinla).*Pnm(:,n-mu+1);
                                        Vfifi=Vfifi+(Cmcosla+Smsinla).*ddPnm(:,n-mu+1);
                                        Vll=Vll+(Cmcosla+Smsinla).*m.^2.*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    VrrH=VrrH.*u+Vrr;
                                    Vrr=0;
                                    VfifiH=VfifiH.*u+Vfifi;
                                    Vfifi=0;
                                    VllH=VllH.*u+Vll;
                                    Vll=0;
                                end
                                
                            elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                                
                                if m<2   
                                    for n=mu:nmax
                                        Smcosla=Sm(n-mu+1).*cosla;
                                        Cmsinla=Cm(n-mu+1).*sinla;

                                        Vrfi=Vrfi+(n+1).*(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*dPnm(:,n-m+1);                                    
                                        Vrl=Vrl+(n+1).*(Smcosla-Cmsinla).*(m*Pnm(:,n-m+1));
                                        Vfil=Vfil+(Smcosla-Cmsinla).*(m*dPnm(:,n-m+1));
                                    end
                                else
                                    for n=mu:nmax
                                        Smcosla=Sm(n-mu+1).*cosla;
                                        Cmsinla=Cm(n-mu+1).*sinla;

                                        Vrfi=Vrfi+(n+1).*(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*dPnm(:,n-mu+1);
                                        Vrl=Vrl+(n+1).*(Smcosla-Cmsinla).*(m*Pnm(:,n-mu+1));
                                        Vfil=Vfil+(Smcosla-Cmsinla).*(m*dPnm(:,n-mu+1));
                                    end
                                end

                                if volbaALFs==2
                                    VrfiH=VrfiH.*u+Vrfi;
                                    Vrfi=0;
                                    VrlH=VrlH.*u+Vrl;
                                    Vrl=0;
                                    VfilH=VfilH.*u+Vfil;
                                    Vfil=0;
                                end
                                
                            elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                                
                                %Coefficients of the ALFs
                                if m==0
                                    anm=sqrt(2)/4.*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)-(m+2)+2);
                                    bnm=((2:nmax)+m+1).*((2:nmax)+m+2)./2./(m+1);
                                    cnm=zeros(nmax-1,1);
                                elseif m==1
                                    anm=1/4.*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)-(m+2)+2);
                                    bnm=((2:nmax)+m+1).*((2:nmax)+m+2)./2./(m+1);
                                    cnm=zeros(nmax-1,1);
                                elseif m==2
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    cnm=sqrt(2)/4.*sqrt((2:nmax).^2-(m-2+1).^2).*sqrt((2:nmax)-(m-2)).*sqrt((2:nmax)+m-2+2);
                                else
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    cnm=1/4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                end
                                  
                                %Spherical harmonic coefficients of orders m+2 a m-2
                                if radenie==0
                                    if m<2
                                        Cplus2=C(index+m+2);
                                        Splus2=S(index+m+2); 
                                        
                                        Cminus2=zeros(nmax-1,1);
                                        Sminus2=zeros(nmax-1,1); 
                                    elseif m>nmax-2
                                        Cplus2=zeros(nmax-m+1,1);
                                        Splus2=zeros(nmax-m+1,1);
                                        
                                        Cminus2=C(index((m-1):end)+m-2);
                                        Sminus2=S(index((m-1):end)+m-2); 
                                    else
                                        Cplus2=C(index((m-1):end)+m+2);
                                        Splus2=S(index((m-1):end)+m+2);
                                        
                                        Cminus2=C(index((m-1):end)+m-2);
                                        Sminus2=S(index((m-1):end)+m-2); 
                                    end
                                elseif radenie==1
                                    if m<2
                                        Cplus2=[zeros(m,1);C((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))];
                                        Splus2=[zeros(m,1);S((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))]; 
                                        
                                        Cminus2=zeros(nmax-1,1);
                                        Sminus2=zeros(nmax-1,1); 
                                    elseif m>nmax-2
                                        if m==2
                                            Cminus2=C(1:(nmax-1));
                                            Sminus2=C(1:(nmax-1))*0;
                                        else
                                            Cminus2=C((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            Sminus2=S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                        
                                        Cplus2=zeros(nmax-m+1,1);
                                        Splus2=zeros(nmax-m+1,1);
                                    else
                                        if m==2
                                            Cminus2=C(1:(nmax-1));
                                            Sminus2=C(1:(nmax-1))*0;
                                        else
                                            Cminus2=C((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            Sminus2=S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                        
                                        Cplus2=[0;0;C((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                        Splus2=[0;0;S((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                    end
                                end
                                
                                %Summation over n                            
                                if m<2   
                                    for n=mu:nmax
                                        %Vzz
                                        Vzz=Vzz+(n+1).*(n+2).*(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);    

                                        %Vxx    
                                        Vxx1=(Cplus2(n-mu+1).*coslaplus2+Splus2(n-mu+1).*sinlaplus2).*anm(n-mu+1);
                                        Vxx2_vnutro=Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla;
                                        Vxx2=Vxx2_vnutro.*(bnm(n-mu+1)-(n+1)*(n+2));
                                        Vxx3=(Cminus2(n-mu+1).*coslaminus2+Sminus2(n-mu+1).*sinlaminus2).*cnm(n-mu+1);
                                        Vxx=Vxx+(Vxx1+Vxx2+Vxx3).*Pnm(:,n-m+1);  

                                        %Vyy
                                        Vyy2=Vxx2_vnutro.*bnm(n-mu+1);
                                        Vyy=Vyy+(Vxx1+Vyy2+Vxx3).*Pnm(:,n-m+1); 
                                    end
                                else
                                    for n=mu:nmax
                                        %Vzz
                                        Vzz=Vzz+(n+1).*(n+2).*(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);

                                        %Vxx
                                        Vxx1=(Cplus2(n-mu+1).*coslaplus2+Splus2(n-mu+1).*sinlaplus2).*anm(n-mu+1);
                                        Vxx2_vnutro=Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla;
                                        Vxx2=Vxx2_vnutro.*(bnm(n-mu+1)-(n+1)*(n+2));
                                        Vxx3=(Cminus2(n-mu+1).*coslaminus2+Sminus2(n-mu+1).*sinlaminus2).*cnm(n-mu+1);
                                        Vxx=Vxx+(Vxx1+Vxx2+Vxx3).*Pnm(:,n-mu+1);

                                        %Vyy
                                        Vyy2=Vxx2_vnutro.*bnm(n-mu+1);
                                        Vyy=Vyy+(Vxx1+Vyy2+Vxx3).*Pnm(:,n-mu+1); 
                                    end
                                end

                                if volbaALFs==2
                                    VzzH=VzzH.*u+Vzz;
                                    Vzz=0;
                                    VxxH=VxxH.*u+Vxx;
                                    Vxx=0;
                                    VyyH=VyyH.*u+Vyy;
                                    Vyy=0;
                                end
                                
                            elseif volbapar(i)==15 %Disturbing tensor Txy_Txz_Tyz
                                
                                %Coefficietns of the ALFs
                                if m==0
                                    dnm=1/4.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt(2).*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)+m+2-2);
                                    gnm=zeros(nmax-1,1);
                                    hnm=gnm;
                                    
                                    betanm=((2:nmax)+2)./2.*sqrt(1+ones(1,nmax-1)).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)-(m+1)+1);
                                    gamanm=gnm;
                                    
                                    minm=((2:nmax)+2)./2.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt(2).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)+m+1-1);
                                    ninm=gnm;
                                elseif m==1
                                    dnm=1/4.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax).^2-(m+2-1).^2).*sqrt((2:nmax)+m+2).*sqrt((2:nmax)+m+2-2);
                                    gnm=-1/4*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax)+1).*sqrt((2:nmax)-1).*((2:nmax)+2);
                                    hnm=zeros(nmax-1,1);
                                    
                                    betanm=((2:nmax)+2)./2.*sqrt((2:nmax)+m+1).*sqrt((2:nmax)-(m+1)+1);
                                    gamanm=-((2:nmax)+2).*sqrt((2:nmax).*((2:nmax)+1)./2);
                                    
                                    minm=((2:nmax)+2)./2.*sqrt((2*(2:nmax)+1)./(2*(2:nmax)-1)).*sqrt((2:nmax)+m+1).*sqrt((2:nmax)+m+1-1);
                                    ninm=gnm;
                                elseif m==2
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=zeros(nmax-1,1);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                elseif m==3
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-3).*sqrt((m:nmax)-2).*sqrt((m:nmax)-1).*sqrt((m:nmax)+2);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                else
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)-(m-2)-2);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                end
                                
                                %Spherical harmonic coefficients of orders m+2, m-2, m+1 and m-1
                                if radenie==0
                                    if m<2
                                        Cplus2=C(index+m+2);
                                        Splus2=S(index+m+2); 
                                        
                                        Cminus2=zeros(nmax-1,1);
                                        Sminus2=zeros(nmax-1,1); 
                                        
                                        Cplus1=C(index+m+1);
                                        Splus1=S(index+m+1);
                                        
                                        if m==0
                                            Cminus1=zeros(nmax-1,1);
                                            Sminus1=zeros(nmax-1,1);
                                        else
                                            Cminus1=C(index+m-1);
                                            Sminus1=S(index+m-1);
                                        end
                                    elseif m>nmax-2
                                        Cplus2=zeros(nmax-m+1,1);
                                        Splus2=zeros(nmax-m+1,1);
                                        
                                        Cminus2=C(index((m-1):end)+m-2);
                                        Sminus2=S(index((m-1):end)+m-2); 
                                        
                                        if m==nmax
                                            Cplus1=zeros(nmax-m+1,1);
                                            Splus1=zeros(nmax-m+1,1);
                                        else
                                            Cplus1=C(index((m-1):end)+m+1);
                                            Splus1=S(index((m-1):end)+m+1);
                                        end
                                        
                                        Cminus1=C(index((m-1):end)+m-1);
                                        Sminus1=S(index((m-1):end)+m-1); 
                                    else
                                        Cplus2=C(index((m-1):end)+m+2);
                                        Splus2=S(index((m-1):end)+m+2);
                                        
                                        Cminus2=C(index((m-1):end)+m-2);
                                        Sminus2=S(index((m-1):end)+m-2); 
                                        
                                        Cplus1=C(index((m-1):end)+m+1);
                                        Splus1=S(index((m-1):end)+m+1);
                                        
                                        Cminus1=C(index((m-1):end)+m-1);
                                        Sminus1=S(index((m-1):end)+m-1); 
                                    end
                                elseif radenie==1
                                    if m<2
                                        Cplus2=[zeros(m,1);C((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))];
                                        Splus2=[zeros(m,1);S((zT+(2*nmaxGGM-2)):(kT+(2*nmaxGGM-2-m)))]; 
                                        
                                        Cminus2=zeros(nmax-1,1);
                                        Sminus2=zeros(nmax-1,1); 
                                        
                                        Cplus1=C((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)));
                                        Splus1=S((zT+(nmaxGGM-1)):(kT+(nmaxGGM-1)));
                                        
                                        if m==1
                                            Cminus1=C(1:(nmax-1));
                                            Sminus1=C(1:(nmax-1))*0;
                                        else                                            
                                            Cminus1=zeros(nmax-1,1);
                                            Sminus1=zeros(nmax-1,1);
                                        end
                                    elseif m>nmax-2
                                        if m==2
                                            Cminus2=C(1:(nmax-1));
                                            Sminus2=C(1:(nmax-1))*0;
                                        else
                                            Cminus2=C((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            Sminus2=S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                        
                                        Cplus2=zeros(nmax-m+1,1);
                                        Splus2=zeros(nmax-m+1,1);
                                        
                                        if m==nmax
                                            Cplus1=zeros(nmax-m+1,1);
                                            Splus1=zeros(nmax-m+1,1);
                                        else
                                            Cplus1=[0;C((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                            Splus1=[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                        end
                                        
                                        Cminus1=C((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                        Sminus1=S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1))); 
                                    else
                                        if m==2
                                            Cminus2=C(1:(nmax-1));
                                            Sminus2=C(1:(nmax-1))*0;
                                        else
                                            Cminus2=C((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                            Sminus2=S((zT-((nmaxGGM-m)*2+3)):(kT-((nmaxGGM-m)*2+3)));
                                        end
                                        
                                        Cplus2=[0;0;C((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                        Splus2=[0;0;S((zT+((nmaxGGM-m)*2+1)):(kT+((nmaxGGM-m)*2-1)))];
                                        
                                        Cplus1=[0;C((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                        Splus1=[0;S((zT+((nmaxGGM-m)+1)):(kT+((nmaxGGM-m))))];
                                        
                                        Cminus1=C((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                        Sminus1=S((zT-((nmaxGGM-m)+1)):(kT-((nmaxGGM-m)+1)));
                                    end
                                end
                                
                                PTnm=[zeros(length(fi),1) Pnm];
                                
                                %Summation over n                            
                                if m<2  
                                    for n=mu:nmax
                                        %Vxy    
                                        Vxy1=(-Cplus2(n-mu+1).*sinlaplus2+Splus2(n-mu+1).*coslaplus2).*dnm(n-mu+1);
                                        Vxy2=(-Cm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*gnm(n-mu+1);
                                        Vxy3=(-Cminus2(n-mu+1).*sinlaminus2+Sminus2(n-mu+1).*coslaminus2).*hnm(n-mu+1);
                                        Vxy=Vxy+(Vxy1+Vxy2+Vxy3).*PTnm(:,n-m+1).*q;  

                                        %Vxz
                                        Vxz1=(Cplus1(n-mu+1).*coslaplus1+Splus1(n-mu+1).*sinlaplus1).*betanm(n-mu+1);
                                        Vxz2=(Cminus1(n-mu+1).*coslaminus1+Sminus1(n-mu+1).*sinlaminus1).*gamanm(n-mu+1);
                                        Vxz=Vxz+(Vxz1+Vxz2).*Pnm(:,n-m+1);

                                        %Vyz
                                        Vyz1=(-Cplus1(n-mu+1).*sinlaplus1+Splus1(n-mu+1).*coslaplus1).*minm(n-mu+1);
                                        Vyz2=(-Cminus1(n-mu+1).*sinlaminus1+Sminus1(n-mu+1).*coslaminus1).*ninm(n-mu+1);
                                        Vyz=Vyz+(Vyz1+Vyz2).*PTnm(:,n-m+1).*q;  
                                    end
                                else
                                    for n=mu:nmax
                                        %Vxy
                                        Vxy1=(-Cplus2(n-mu+1).*sinlaplus2+Splus2(n-mu+1).*coslaplus2).*dnm(n-mu+1);
                                        Vxy2=(-Cm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*gnm(n-mu+1);
                                        Vxy3=(-Cminus2(n-mu+1).*sinlaminus2+Sminus2(n-mu+1).*coslaminus2).*hnm(n-mu+1);
                                        Vxy=Vxy+(Vxy1+Vxy2+Vxy3).*PTnm(:,n-mu+1).*q;

                                        %Vxz
                                        Vxz1=(Cplus1(n-mu+1).*coslaplus1+Splus1(n-mu+1).*sinlaplus1).*betanm(n-mu+1);
                                        Vxz2=(Cminus1(n-mu+1).*coslaminus1+Sminus1(n-mu+1).*sinlaminus1).*gamanm(n-mu+1);
                                        Vxz=Vxz+(Vxz1+Vxz2).*Pnm(:,n-mu+1);

                                        %Vyz
                                        Vyz1=(-Cplus1(n-mu+1).*sinlaplus1+Splus1(n-mu+1).*coslaplus1).*minm(n-mu+1);
                                        Vyz2=(-Cminus1(n-mu+1).*sinlaminus1+Sminus1(n-mu+1).*coslaminus1).*ninm(n-mu+1);
                                        Vyz=Vyz+(Vyz1+Vyz2).*PTnm(:,n-mu+1).*q;  
                                    end
                                end

                                if volbaALFs==2
                                    VxyH=VxyH.*u+Vxy;
                                    Vxy=0;
                                    VxzH=VxzH.*u+Vxz;
                                    Vxz=0;
                                    VyzH=VyzH.*u+Vyz;
                                    Vyz=0;
                                end
                                
                            elseif volbapar(i)==16 %Gravity
                               
                                if m<2  
                                    for n=mu:nmax 
                                        Cmcosla=Cm(n-mu+1).*cosla;
                                        Smsinla=Sm(n-mu+1).*sinla;
                                        
                                        Wr=Wr+(n+1).*(Cmcosla+Smsinla).*Pnm(:,n-m+1); 
                                        Wlambda=Wlambda+(-Cm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*(m*Pnm(:,n-m+1)); 
                                        Wfi=Wfi+(Cmcosla+Smsinla).*dPnm(:,n-m+1); 
                                    end
                                else
                                    for n=mu:nmax 
                                        Cmcosla=Cm(n-mu+1).*cosla;
                                        Smsinla=Sm(n-mu+1).*sinla;
                                        
                                        Wr=Wr+(n+1).*(Cmcosla+Smsinla).*Pnm(:,n-mu+1);
                                        Wlambda=Wlambda+(-Cm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*(m*Pnm(:,n-mu+1)); 
                                        Wfi=Wfi+(Cmcosla+Smsinla).*dPnm(:,n-mu+1); 
                                    end
                                end

                                if volbaALFs==2
                                    WrH=WrH.*u+Wr;
                                    Wr=0;
                                    WlambdaH=WlambdaH.*u+Wlambda;
                                    Wlambda=0;
                                    WfiH=WfiH.*u+Wfi;
                                    Wfi=0;
                                end

                            elseif volbapar(i)==17 %Gravity sa
                            
                                if m<2    
                                    for n=mu:nmax  
                                        g_sa=g_sa+(n+1).*(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);  
                                    end
                                else
                                    for n=mu:nmax  
                                        g_sa=g_sa+(n+1).*(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    g_saH=g_saH.*u+g_sa;
                                    g_sa=0;
                                end

                            elseif volbapar(i)==18 %Gravity potential
                                                               
                                if m<2 
                                    for n=mu:nmax 
                                        W=W+(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);    
                                    end
                                else
                                    for n=mu:nmax 
                                        W=W+(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    WH=WH.*u+W;
                                    W=0;
                                end

                            elseif volbapar(i)==19 %Gravity anomaly sa  
                                                               
                                if m<2    
                                    for n=mu:nmax 
                                        anomalia_sa=anomalia_sa+(n-1).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);    
                                    end
                                else
                                    for n=mu:nmax 
                                        anomalia_sa=anomalia_sa+(n-1).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    anomalia_saH=anomalia_saH.*u+anomalia_sa;
                                    anomalia_sa=0;
                                end

                            elseif volbapar(i)==20 %Gravity disturbance
                              
                                if m<2    
                                    for n=mu:nmax 
                                        Cmcosla=Cm(n-mu+1).*cosla;
                                        Smsinla=Sm(n-mu+1).*sinla;
                                        
                                        Wr=Wr+(n+1).*(Cmcosla+Smsinla).*Pnm(:,n-m+1); 
                                        Wlambda=Wlambda+(-Cm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*(m*Pnm(:,n-m+1)); 
                                        Wfi=Wfi+(Cmcosla+Smsinla).*dPnm(:,n-m+1); 
                                    end
                                else
                                    for n=mu:nmax 
                                        Cmcosla=Cm(n-mu+1).*cosla;
                                        Smsinla=Sm(n-mu+1).*sinla;
                                        
                                        Wr=Wr+(n+1).*(Cmcosla+Smsinla).*Pnm(:,n-mu+1);
                                        Wlambda=Wlambda+(-Cm(n-mu+1).*sinla+Sm(n-mu+1).*cosla).*(m*Pnm(:,n-mu+1)); 
                                        Wfi=Wfi+(Cmcosla+Smsinla).*dPnm(:,n-mu+1);  
                                    end
                                end

                                if m==0
                                    for n=mu:10
                                        Ur=Ur+(n+1).*CElm(n-mu+1).*cosla.*Pnm(:,n-m+1); 
                                        Ufi=Ufi+CElm(n-mu+1).*cosla.*dPnm(:,n-m+1); 
                                    end
                                end

                                if volbaALFs==2
                                    WrH=WrH.*u+Wr;
                                    Wr=0;
                                    WlambdaH=WlambdaH.*u+Wlambda;
                                    Wlambda=0;
                                    WfiH=WfiH.*u+Wfi;
                                    Wfi=0;
                                    UrH=UrH.*u+Ur;
                                    Ur=0;
                                    UfiH=UfiH.*u+Ufi;
                                    Ufi=0;
                                end

                            elseif volbapar(i)==21 %Gravity disturbance sa
                             
                                if m<2     
                                    for n=mu:nmax 
                                        porucha_sa=porucha_sa+(n+1).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1); 
                                    end
                                else
                                    for n=mu:nmax 
                                        porucha_sa=porucha_sa+(n+1).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    porucha_saH=porucha_saH.*u+porucha_sa;
                                    porucha_sa=0;
                                end

                            elseif volbapar(i)==22 %Height anomaly Ell
                                                             
                                if m<2        
                                    for n=mu:nmax 
                                        zetaEl=zetaEl+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);                                    
                                    end
                                else
                                    for n=mu:nmax 
                                        zetaEl=zetaEl+(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    zetaElH=zetaElH.*u+zetaEl;
                                    zetaEl=0;
                                end
                                
                            elseif volbapar(i)==23 %Height anomaly
                                
                                if m<2 %Spherical harmonics of geoid undulation of orders 0 a 1
                                    for nn=m:1
                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFFIED
                                        % fnALFS
                                        
                                        zeta_H=zeta_H+1./(R./r).^nn.*(HCm(nn-m+1).*cosla+HSm(nn-m+1).*sinla).*Pnm(:,nn-m+1); 
                                    end
                                end
                                
                                                               
                                if m<2   
                                    for n=mu:nmax 
                                        deltaCmcosla_Smsinla_Pnm=(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);
                                        
                                        zeta_N1c=zeta_N1c+deltaCmcosla_Smsinla_Pnm;  
                                        
                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFFIED
                                        % fnALFS
                                        
                                        zeta_H=zeta_H+1./(R./r).^n.*(HCm(n-mu+3-m).*cosla+HSm(n-mu+3-m).*sinla).*Pnm(:,n-m+1); 
                                        zeta_dg=zeta_dg+(n+1).*deltaCmcosla_Smsinla_Pnm;  
                                    end
                                else
                                    for n=mu:nmax 
                                        deltaCmcosla_Smsinla_Pnm=(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                        
                                        zeta_N1c=zeta_N1c+deltaCmcosla_Smsinla_Pnm;
                                        
                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFFIED
                                        % fnALFS
                                        
                                        zeta_H=zeta_H+1./(R./r).^n.*(HCm(n-mu+1).*cosla+HSm(n-mu+1).*sinla).*Pnm(:,n-mu+1); 
                                        zeta_dg=zeta_dg+(n+1).*deltaCmcosla_Smsinla_Pnm;
                                    end
                                end
                                    
                                if volbaALFs==2
                                    zeta_N1cH=zeta_N1cH.*u+zeta_N1c;
                                    zeta_N1c=0;
                                    zeta_HH=zeta_HH.*u+zeta_H;
                                    zeta_H=0;
                                    zeta_dgH=zeta_dgH.*u+zeta_dg;
                                    zeta_dg=0;
                                end
                                
                            elseif volbapar(i)==24 %Second radial derivative of disturbing potential
                               
                                if m<2     
                                    for n=mu:nmax    
                                        T_rr=T_rr+(n+1).*(n+2).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);    
                                    end
                                else
                                    for n=mu:nmax    
                                        T_rr=T_rr+(n+1).*(n+2).*(deltaCm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                    end
                                end

                                if volbaALFs==2
                                    T_rrH=T_rrH.*u+T_rr;
                                    T_rr=0;
                                end

                            elseif volbapar(i)==25 %Second radial derivative of gravity potential
                                                              
                               if m<2
                                   for n=mu:nmax 
                                        Wrr=Wrr+(n+1).*(n+2).*(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-m+1);     
                                   end
                               else
                                   for n=mu:nmax 
                                        Wrr=Wrr+(n+1).*(n+2).*(Cm(n-mu+1).*cosla+Sm(n-mu+1).*sinla).*Pnm(:,n-mu+1);
                                   end
                               end

                               if volbaALFs==2
                                   WrrH=WrrH.*u+Wrr;
                                   Wrr=0;
                               end
                                    
                            end                                              
                        end


                    end
                    
                    %Update of the progress bar
                    set(progressbar,'string',...
                        'Progress: Matrix multiplications...',...
                        'fontsize',8); drawnow;

                    clear Pnm dPnm ddPnm Cm Sm C CElm deltaC deltaCm ...
                        S u t q q2 index sinla cosla ...
                        deltaCmcosla_Smsinla_Pnm Cmcosla deltaCmcosla ...
                        Smcosla deltaCmsinla Smsinla enm singdPnm ...
                        qu tu singddPnm 

                    %Computation of the normal gravity for eta, xi, Theta,
                    %N, zeta el, zeta
                    if any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==10) || any(volbapar==22) || any(volbapar==23)                        
                        bEl=aEl*sqrt(1-eEl^2);
                        EEl=sqrt(aEl^2-bEl^2);
                        
                        %Computation of ellipsoidal harmonic coordinates
                        [X,Y,Z]=geodetic2ecef(fi,lambda,h,[aEl eEl]');
                        ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                        betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                        clear X Y Z
                        
                        wgama=sqrt((ugama.^2+EEl^2*sin(betagama).^2)./(ugama.^2+EEl^2));
                        qgama=1/2*((1+(3*ugama.^2)./EEl^2).*atan(EEl./ugama)-3*ugama./EEl);
                        qgama_=3*(1+(ugama.^2)./EEl^2).*(1-ugama./EEl.*atan(EEl./ugama))-1;
                        qgama0=1/2*((1+(3*bEl^2)/EEl^2)*atan(EEl/bEl)-3*bEl/EEl);
                        
                        gamau=-1./wgama.*(GMEl./(ugama.^2+EEl.^2)+(omegaEl.^2.*aEl.^2.*EEl)./(ugama.^2+EEl.^2).*(qgama_./qgama0).*(1./2.*sin(betagama).^2-1./6)-omegaEl.^2.*ugama.*cos(betagama).^2);
                        gamabeta=-1./wgama.*(-(omegaEl.^2.*aEl.^2.*qgama)./(sqrt(ugama.^2+EEl.^2).*qgama0)+omegaEl.^2.*sqrt(ugama.^2+EEl.^2)).*sin(betagama).*cos(betagama);
                        
                        clear ugama betagama wgama qgama qgama_ qgama0
                        
                        gamaP=sqrt(gamau.^2+gamabeta.^2);
                        
                        clear gamau gamabeta
                    end
                    
                    %% Final computation of functionals of the geopotential 
                    for i=1:pocetpar 
                        if volbapar(i)==1                
                        elseif volbapar(i)==2 %Deflection of the vertical eta

                            if volbaALFs==1 || volbaALFs==3
                                Pg=-(GM./(r.^2.*gamaP.*cos(fiG)).*(sum(eta,2)))*(180/pi)*3600;
                            elseif volbaALFs==2
                                Pg=-(GM./(r.^2.*gamaP.*cos(fiG)).*(etaH*1e280))*(180/pi)*3600;
                            end

                            Pg(fi==pi/2 | fi==-pi/2)=0;
                            
                            clear eta etaH

                        elseif volbapar(i)==3 %Deflection of the vertical xi

                            if volbaALFs==1 || volbaALFs==3
                                Pg=-(GM./(r.^2.*gamaP).*(sum(ksi,2)))*(180/pi)*3600;
                            elseif volbaALFs==2
                                Pg=-(GM./(r.^2.*gamaP).*(ksiH*1e280))*(180/pi)*3600;
                            end

                            clear ksi ksiH

                        elseif volbapar(i)==4 %Deflection of the vertical Theta

                            if volbaALFs==1 || volbaALFs==3
                                Teta=-(GM./(r.^2.*gamaP.*cos(fiG)).*(sum(Teta,2)))*(180/pi)*3600;
                                Tksi=-(GM./(r.^2.*gamaP).*(sum(Tksi,2)))*(180/pi)*3600;
                                
                                Talfa=atan2(Teta,Tksi);
                                Talfa(Talfa<0)=Talfa(Talfa<0)+2*pi;
                            
                                Pg=[sqrt(Teta.^2+Tksi.^2) rad2deg(Talfa)];
                            
                                clear Teta Tksi Talfa
                                
                                Pg(fi==pi/2 | fi==-pi/2,:)=0;
                                
                            elseif volbaALFs==2
                                TetaH=-(GM./(r.^2.*gamaP.*cos(fiG)).*(TetaH*1e280))*(180/pi)*3600;
                                TksiH=-(GM./(r.^2.*gamaP).*(TksiH*1e280))*(180/pi)*3600;
                                
                                Talfa=atan2(TetaH,TksiH);
                                Talfa(Talfa<0)=Talfa(Talfa<0)+2*pi;

                                Pg=[sqrt(TetaH.^2+TksiH.^2) rad2deg(Talfa)];

                                clear TetaH TksiH Talfa
                                
                                Pg(fi==pi/2 | fi==-pi/2,:)=0;
                            end

                        elseif volbapar(i)==5 %Disturbing potential

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.*(sum(T,2)));
                            elseif volbaALFs==2
                                Pg=(GM./r.*(TH*1e280));
                            end

                            clear T TH

                        elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                            
                            if volbaALFs==1 || volbaALFs==3
                                Trr=(GM./r.^3.*(sum(Trr,2)))*10^9;
                                Tfifi=(GM./r.^3.*(sum(Tfifi,2)))*10^9;
                                Tll=(GM./(r.^3.*cos(fiG).^2).*(sum(Tll,2)))*10^9;
                                
                                Pg=Trr;
                                clear Trr
                                Pg=[Pg Tfifi];
                                clear Tfifi
                                Pg=[Pg -Tll];
                                clear Tll
                            elseif volbaALFs==2
                                TrrH=(GM./r.^3.*(TrrH*1e280))*10^9;
                                TfifiH=(GM./r.^3.*(TfifiH*1e280))*10^9;
                                TllH=(GM./(r.^3.*cos(fiG).^2).*(TllH*1e280))*10^9;
                                
                                Pg=TrrH;
                                clear TrrH
                                Pg=[Pg TfifiH];
                                clear TfifiH
                                Pg=[Pg -TllH];
                                clear TllH
                            end
                            
                        elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                            
                            if volbaALFs==1 || volbaALFs==3
                                Trfi=(GM./r.^3.*(sum(Trfi,2)))*10^9;
                                Trl=(GM./(r.^3.*cos(fiG)).*(sum(Trl,2)))*10^9;
                                Tfil=(GM./(r.^3.*cos(fiG)).*(sum(Tfil,2)))*10^9;
                                
                                Pg=-Trfi;
                                clear Trfi
                                Pg=[Pg -Trl];
                                clear Trl
                                Pg=[Pg Tfil];
                                clear Tfil
                            elseif volbaALFs==2
                                TrfiH=(GM./r.^3.*(TrfiH*1e280))*10^9;
                                TrlH=(GM./(r.^3.*cos(fiG)).*(TrlH*1e280))*10^9;
                                TfilH=(GM./(r.^3.*cos(fiG)).*(TfilH*1e280))*10^9;
                                
                                Pg=-TrfiH;
                                clear TrfiH
                                Pg=[Pg -TrlH];
                                clear TrlH
                                Pg=[Pg TfilH];
                                clear TfilH
                            end
                            
                        elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                            
                            clear deltaCplus2 Splus2 coslaplus2 sinlaplus2 ...
                                deltaCminus2 Sminus2 coslaminus2 sinlaminus2 ...
                                Txx1 Txx2 Txx3 Tyy2 Txx2_vnutro anm bnm cnm
                            
                            if volbaALFs==1 || volbaALFs==3
                                Tzz=(GM./r.^3.*(sum(Tzz,2)))*10^9;
                                Txx=(GM./r.^3.*(sum(Txx,2)))*10^9;
                                Tyy=(GM./r.^3.*(sum(Tyy,2)))*10^9;
                                
                                Pg=[Txx];
                                clear Txx
                                Pg=[Pg -Tyy];
                                clear Tyy
                                Pg=[Pg Tzz];
                                clear Tzz
                            elseif volbaALFs==2
                                TzzH=(GM./r.^3.*(TzzH*1e280))*10^9;
                                TxxH=(GM./r.^3.*(TxxH*1e280))*10^9;
                                TyyH=(GM./r.^3.*(TyyH*1e280))*10^9;
                                
                                Pg=[TxxH];
                                clear TxxH
                                Pg=[Pg -TyyH];
                                clear TyyH
                                Pg=[Pg TzzH];
                                clear TzzH
                            end
                            
                        elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                            
                            clear deltaCplus2 Splus2 coslaplus2 sinlaplus2 ...
                                deltaCminus2 Sminus2 coslaminus2 sinlaminus2 ...
                                deltaCplus1 Splus1 coslaplus1 sinlaplus1 ...
                                deltaCminus1 Sminus1 coslaminus1 sinlaminus1 ...
                                PTnm Txy1 Txy2 Txy3 Txz1 Txz2 Tyz1 Tyz2 ...
                                dnm gnm hnm betanm gamanm minm ninm
                            
                            if volbaALFs==1 || volbaALFs==3
                                Txy=(GM./r.^3.*(sum(Txy,2)))*10^9;
                                Txz=(GM./r.^3.*(sum(Txz,2)))*10^9;
                                Tyz=(GM./r.^3.*(sum(Tyz,2)))*10^9;
                                
                                Pg=[Txy];
                                clear Txy
                                Pg=[Pg Txz];
                                clear Txz
                                Pg=[Pg Tyz];
                                clear Tyz
                            elseif volbaALFs==2
                                TxyH=(GM./r.^3.*(TxyH*1e280))*10^9;
                                TxzH=(GM./r.^3.*(TxzH*1e280))*10^9;
                                TyzH=(GM./r.^3.*(TyzH*1e280))*10^9;
                                
                                Pg=[TxyH];
                                clear TxyH
                                Pg=[Pg TxzH];
                                clear TxzH
                                Pg=[Pg TyzH];
                                clear TyzH
                            end

                        elseif volbapar(i)==10 %Geoid undulation
                            
                            clear HC HCm HS HSm
                            
                            G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                            ro=2670; %Density of the crust

                            if volbaALFs==1 || volbaALFs==3
                                N1c=(GM./r.*(sum(N1c,2)))./gamaP;
                                H(H<0)=H(H<0)*0; %H is set to zero in the areas of oceans and seas
                                Pg=N1c-(2*pi*G*ro*H.^2)./gamaP;
                                
                                clear N1c H
                                
                            elseif volbaALFs==2
                                N1cH=(GM./r.*(N1cH*1e280))./gamaP;
                                HH=HH*1e280;
                                HH(HH<0)=HH(HH<0)*0;
                                Pg=N1cH-(2*pi*G*ro*HH.^2)./gamaP;
                                
                                clear N1cH HH
                            end
                             
                        elseif volbapar(i)==11 %Gravitational potential

                            if volbaALFs==1 || volbaALFs==3
                                if nmin0==1 %Zero degree term
                                    Pg=(GM./r.*(1+sum(V,2)));
                                elseif nmin0==0
                                    Pg=(GM./r.*(sum(V,2)));
                                end
                            elseif volbaALFs==2
                                if nmin0==1 %Zero degree term
                                    Pg=(GM./r.*(1+VH*1e280));
                                elseif nmin0==0
                                    Pg=(GM./r.*(VH*1e280));
                                end
                            end

                            clear V VH

                        elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                            
                            if volbaALFs==1 || volbaALFs==3
                                if nmin0==1 %Zero degree term
                                    Vrr=(GM./r.^3.*(2+sum(Vrr,2)))*10^9;
                                elseif nmin0==0
                                    Vrr=(GM./r.^3.*(sum(Vrr,2)))*10^9;
                                end
                                Vfifi=(GM./r.^3.*(sum(Vfifi,2)))*10^9;
                                Vll=(GM./(r.^3.*cos(fiG).^2).*(sum(Vll,2)))*10^9;
                                
                                Pg=Vrr;
                                clear Vrr
                                Pg=[Pg Vfifi];
                                clear Vfifi
                                Pg=[Pg -Vll];
                                clear Vll
                            elseif volbaALFs==2
                                if nmin0==1 %Zero degree term
                                    VrrH=(GM./r.^3.*(2+VrrH*1e280))*10^9;
                                elseif nmin0==0
                                    VrrH=(GM./r.^3.*(VrrH*1e280))*10^9;
                                end
                                VfifiH=(GM./r.^3.*(VfifiH*1e280))*10^9;
                                VllH=(GM./(r.^3.*cos(fiG).^2).*(VllH*1e280))*10^9;
                                
                                Pg=VrrH;
                                clear VrrH
                                Pg=[Pg VfifiH];
                                clear VfifiH
                                Pg=[Pg -VllH];
                                clear VllH
                            end
                            
                        elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vrfi=(GM./r.^3.*(sum(Vrfi,2)))*10^9;
                                Vrl=(GM./(r.^3.*cos(fiG)).*(sum(Vrl,2)))*10^9;
                                Vfil=(GM./(r.^3.*cos(fiG)).*(sum(Vfil,2)))*10^9;
                                
                                Pg=-Vrfi;
                                clear Vrfi
                                Pg=[Pg -Vrl];
                                clear Vrl
                                Pg=[Pg Vfil];
                                clear Vfil
                            elseif volbaALFs==2
                                VrfiH=(GM./r.^3.*(VrfiH*1e280))*10^9;
                                VrlH=(GM./(r.^3.*cos(fiG)).*(VrlH*1e280))*10^9;
                                VfilH=(GM./(r.^3.*cos(fiG)).*(VfilH*1e280))*10^9;
                                
                                Pg=-VrfiH;
                                clear VrfiH
                                Pg=[Pg -VrlH];
                                clear VrlH
                                Pg=[Pg VfilH];
                                clear VfilH
                            end
                            
                        elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                            
                            clear Cplus2 Splus2 coslaplus2 sinlaplus2 ...
                                Cminus2 Sminus2 coslaminus2 sinlaminus2 ...
                                Vxx1 Vxx2 Vxx3 Vyy2 Vxx2_vnutro anm bnm cnm
                            
                            if volbaALFs==1 || volbaALFs==3
                                if nmin0==1 %Zero degree term
                                    Vzz=(GM./r.^3.*(2+sum(Vzz,2)))*10^9;
                                elseif nmin0==0
                                    Vzz=(GM./r.^3.*(sum(Vzz,2)))*10^9;
                                end
                                
                                if nmin0==1 %Zero degree term
                                    Vxx=(GM./r.^3.*(-1+sum(Vxx,2)))*10^9;
                                elseif nmin0==0
                                    Vxx=(GM./r.^3.*(sum(Vxx,2)))*10^9;
                                end
                                
                                if nmin0==1 %Zero degree term
                                    Vyy=(GM./r.^3.*(1+sum(Vyy,2)))*10^9;
                                elseif nmin0==0
                                    Vyy=(GM./r.^3.*(sum(Vyy,2)))*10^9;
                                end
                                
                                Pg=[Vxx];
                                clear Vxx
                                Pg=[Pg -Vyy];
                                clear Vyy
                                Pg=[Pg Vzz];
                                clear Vzz
                            elseif volbaALFs==2
                                if nmin0==1 %Zero degree term
                                    VzzH=(GM./r.^3.*(2+VzzH*1e280))*10^9;
                                elseif nmin0==0
                                    VzzH=(GM./r.^3.*(VzzH*1e280))*10^9;
                                end
                                
                                if nmin0==1 %Zero degree term
                                    VxxH=(GM./r.^3.*(-1+VxxH*1e280))*10^9;
                                elseif nmin0==0
                                    VxxH=(GM./r.^3.*(VxxH*1e280))*10^9;
                                end
                                
                                if nmin0==1 %Zero degree term
                                    VyyH=(GM./r.^3.*(1+VyyH*1e280))*10^9;
                                elseif nmin0==0
                                    VyyH=(GM./r.^3.*(VyyH*1e280))*10^9;
                                end
                                
                                Pg=[VxxH];
                                clear VxxH
                                Pg=[Pg -VyyH];
                                clear VyyH
                                Pg=[Pg VzzH];
                                clear VzzH
                            end
                            
                        elseif volbapar(i)==15 %Vxy_Vxz_Vyz
                            
                            clear Cplus2 Splus2 coslaplus2 sinlaplus2 ...
                                Cminus2 Sminus2 coslaminus2 sinlaminus2 ...
                                Cplus1 Splus1 coslaplus1 sinlaplus1 ...
                                Cminus1 Sminus1 coslaminus1 sinlaminus1 ...
                                PTnm Vxy1 Vxy2 Vxy3 Vxz1 Vxz2 Vyz1 Vyz2 ...
                                dnm gnm hnm betanm gamanm minm ninm
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vxy=(GM./r.^3.*(sum(Vxy,2)))*10^9;
                                Vxz=(GM./r.^3.*(sum(Vxz,2)))*10^9;
                                Vyz=(GM./r.^3.*(sum(Vyz,2)))*10^9;
                                
                                Pg=[Vxy];
                                clear Vxy
                                Pg=[Pg Vxz];
                                clear Vxz
                                Pg=[Pg Vyz];
                                clear Vyz
                            elseif volbaALFs==2
                                VxyH=(GM./r.^3.*(VxyH*1e280))*10^9;
                                VxzH=(GM./r.^3.*(VxzH*1e280))*10^9;
                                VyzH=(GM./r.^3.*(VyzH*1e280))*10^9;
                                
                                Pg=[VxyH];
                                clear VxyH
                                Pg=[Pg VxzH];
                                clear VxzH
                                Pg=[Pg VyzH];
                                clear VyzH
                            end
                            
                        elseif volbapar(i)==16 %Gravity

                            if volbaALFs==1 || volbaALFs==3
                                if nmin0==1 %Zero degree term
                                    Wr=(-GM./r.^2.*(1+sum(Wr,2))+omegaEl^2.*r.*(cos(fiG).^2));
                                elseif nmin0==0
                                    Wr=(-GM./r.^2.*(sum(Wr,2))+omegaEl^2.*r.*(cos(fiG).^2));
                                end
                                Wlambda=GM./r.*Wlambda;
                                Wfi=(GM./r.*(sum(Wfi,2))-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                
                                Pg=sqrt(Wr.^2+(Wlambda./(r.*cos(fiG))).^2+(1./r.*Wfi).^2)*10^5;
                                
                                clear Wr Wlambda Wfi
                            elseif volbaALFs==2
                                if nmin0==1 %Zero degree term
                                    WrH=(-GM./r.^2.*(1+WrH*1e280)+omegaEl^2.*r.*(cos(fiG).^2));
                                elseif nmin0==0
                                    WrH=(-GM./r.^2.*(WrH*1e280)+omegaEl^2.*r.*(cos(fiG).^2));
                                end
                                WlambdaH=GM./r.*WlambdaH*1e280;
                                WfiH=(GM./r.*WfiH*1e280-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                
                                Pg=sqrt(WrH.^2+(WlambdaH./(r.*cos(fiG))).^2+(1./r.*WfiH).^2)*10^5;
                                
                                clear WrH WlambdaH WfiH
                            end

                        elseif volbapar(i)==17 %Gravity sa

                            if volbaALFs==1 || volbaALFs==3
                                if nmin0==1 %Zero degree term
                                    g_sa=(-GM./r.^2.*(1+sum(g_sa,2))+omegaEl^2.*r.*(cos(fiG).^2));
                                elseif nmin0==0
                                    g_sa=(-GM./r.^2.*(sum(g_sa,2))+omegaEl^2.*r.*(cos(fiG).^2));
                                end
                                Pg=sqrt(g_sa.^2)*10^5;
                                
                                clear g_sa
                            elseif volbaALFs==2
                                if nmin0==1 %Zero degree term
                                    g_saH=(-GM./r.^2.*(1+g_saH*1e280)+omegaEl^2.*r.*(cos(fiG).^2));
                                elseif nmin0==0
                                    g_saH=(-GM./r.^2.*(g_saH*1e280)+omegaEl^2.*r.*(cos(fiG).^2));
                                end
                                Pg=sqrt(g_saH.^2)*10^5;
                                
                                clear g_saH
                            end

                        elseif volbapar(i)==18 %Gravity potential

                            if volbaALFs==1 || volbaALFs==3
                                if nmin0==1 %Zero degree term
                                    Pg=(GM./r.*(1+sum(W,2)))+1/2*omegaEl.^2.*r.^2.*cos(fiG).^2;
                                elseif nmin0==0
                                    Pg=(GM./r.*(sum(W,2)))+1/2*omegaEl.^2.*r.^2.*cos(fiG).^2;
                                end
                            elseif volbaALFs==2
                                if nmin0==1 %Zero degree term
                                    Pg=(GM./r.*(1+WH*1e280))+1/2*omegaEl.^2.*r.^2.*cos(fiG).^2;
                                elseif nmin0==0
                                    Pg=(GM./r.*(WH*1e280))+1/2*omegaEl.^2.*r.^2.*cos(fiG).^2;
                                end
                            end

                            clear W WH

                        elseif volbapar(i)==19 %Gravity anomaly sa

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.^2.*(sum(anomalia_sa,2)))*10^5;
                            elseif volbaALFs==2
                                Pg=(GM./r.^2.*(anomalia_saH*1e280))*10^5;
                            end

                            clear anomalia_sa anomalia_saH
                            
                        elseif volbapar(i)==20 %Gravity disturbance

                            if volbaALFs==1 || volbaALFs==3
                                Wr=(-GM./r.^2.*(1+sum(Wr,2))+omegaEl^2.*r.*(cos(fiG).^2));
                                Wlambda=GM./r.*Wlambda;
                                Wfi=(GM./r.*(sum(Wfi,2))-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                Ur=(-GM./r.^2.*(1+sum(Ur,2))+omegaEl^2.*r.*(cos(fiG).^2));
                                Ufi=(GM./r.*(sum(Ufi,2))-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                
                                g=sqrt(Wr.^2+(Wlambda./(r.*cos(fiG))).^2+(1./r.*Wfi).^2);
                                clear Wr Wlambda Wfi
                                
                                gama_GGM=sqrt(Ur.^2+(1./r.*Ufi).^2);
                                clear Ur Ufi

                                Pg=(g-gama_GGM)*10^5;
                                clear g gama_GGM
                            elseif volbaALFs==2
                                WrH=(-GM./r.^2.*(1+WrH*1e280)+omegaEl^2.*r.*(cos(fiG).^2));
                                WlambdaH=GM./r.*WlambdaH*1e280;
                                WfiH=(GM./r.*WfiH*1e280-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                UrH=(-GM./r.^2.*(1+UrH*1e280)+omegaEl^2.*r.*(cos(fiG).^2));
                                UfiH=(GM./r.*UfiH*1e280-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                
                                g=sqrt(WrH.^2+(WlambdaH./(r.*cos(fiG))).^2+(1./r.*WfiH).^2);
                                clear WrH WlambdaH WfiH
                                
                                gama_GGM=sqrt(UrH.^2+(1./r.*UfiH).^2);
                                clear UrH UfiH

                                Pg=(g-gama_GGM)*10^5;
                                clear g gama_GGM
                            end
                            
                        elseif volbapar(i)==21 %Gravity disturbance sa

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.^2.*(sum(porucha_sa,2)))*10^5;
                            elseif volbaALFs==2
                                Pg=(GM./r.^2.*(porucha_saH*1e280))*10^5;
                            end

                            clear porucha_sa porucha_saH
                        elseif volbapar(i)==22 %Height anomaly Ell  

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./(r.*gamaP).*(sum(zetaEl,2)));
                            elseif volbaALFs==2
                                Pg=(GM./(r.*gamaP).*(zetaElH*1e280));
                            end

                            clear zetaEl zetaElH

                        elseif volbapar(i)==23 %Height anomaly
                            
                            clear HC HCm HS HSm
                            
                            G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                            ro=2670; %Density of the crust

                            if volbaALFs==1 || volbaALFs==3
                                zeta_zetaEl=(GM./r.*(sum(zeta_N1c,2)))./gamaP;
                                zeta_N1c=(GM./r.*(sum(zeta_N1c,2)))./gamaP;
                                zeta_dg=(GM./r.^2.*(sum(zeta_dg,2)));
                                
                                zeta_H(zeta_H<0)=zeta_H(zeta_H<0)*0; %H is set to zero in the areas of oceans and seas                           
                                zeta_N=zeta_N1c-(2*pi*G*ro*zeta_H.^2)./gamaP;

                                Pg=zeta_zetaEl-zeta_dg.*((zeta_H+zeta_N)./gamaP);
                            
                                clear zeta_zetaEl zeta_dg zeta_H zeta_N zeta_N1c
                            elseif volbaALFs==2
                                zeta_zetaElH=(GM./r.*(zeta_N1cH*1e280))./gamaP;
                                zeta_N1cH=(GM./r.*(zeta_N1cH*1e280))./gamaP;
                                zeta_dgH=(GM./r.^2.*(zeta_dgH*1e280));
                                zeta_HH=zeta_HH*1e280;
                                
                                zeta_HH(zeta_HH<0)=zeta_HH(zeta_HH<0)*0;                          
                                zeta_NH=zeta_N1cH-(2*pi*G*ro*zeta_HH.^2)./gamaP;

                                Pg=zeta_zetaElH-zeta_dgH.*((zeta_HH+zeta_NH)./gamaP);
                            
                                clear zeta_zetaElH zeta_dgH zeta_HH zeta_NH zeta_N1cH
                            end
                            
                        elseif volbapar(i)==24 %Second radial derivative of disturbing potential

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.^3.*(sum(T_rr,2)))*10^9;
                            elseif volbaALFs==2
                                Pg=(GM./r.^3.*(T_rrH*1e280))*10^9;
                            end

                            clear T_rr T_rrH

                        elseif volbapar(i)==25 %Second radial derivative of gravity potential

                            if volbaALFs==1 || volbaALFs==3
                                if nmin0==1 %Zero degree term
                                    Pg=(GM./r.^3.*(2+sum(Wrr,2))+omegaEl^2.*cos(fiG).^2)*10^9;
                                elseif nmin0==0
                                    Pg=(GM./r.^3.*(sum(Wrr,2))+omegaEl^2.*cos(fiG).^2)*10^9;
                                end
                            elseif volbaALFs==2
                                if nmin0==1 %Zero degree term
                                    Pg=(GM./r.^3.*(2+WrrH*1e280)+omegaEl^2.*cos(fiG).^2)*10^9;
                                elseif nmin0==0
                                    Pg=(GM./r.^3.*(WrrH*1e280)+omegaEl^2.*cos(fiG).^2)*10^9;
                                end
                            end

                            clear Wrr WrrH
                        end

                        if i==1
                            P=Pg;
                            clear Pg
                        else
                            P=[P Pg];
                            if i==pocetpar
                                clear Pg
                            end
                        end
                    end 
                    
                    clear q q2 gamaP
                    
                    %Update of the progress bar
                    set(progressbar,'string','','fontsize',8); drawnow;
                    
                end
            end

            %% Computation of commission error            
            if STD==1
                
                tic %Start clock to measure computation time
                
                %Warning for the computation of modified fnALFs
                if volbaALFs==2 || volbaALFs==3
                    volbaALFs=1;
                    warn1=warndlg('In order to compute the commission error of functionals of the geopotential, only the standard forward column method can be used for the evaluation of fnALFs. After clicking OK, this method will be applied in the following computations.');
                    waitfor(warn1);
                end
                
                %Check if commission error of tensors is to be computed in
                %the LNOF (not allowed)
                if any(volbapar==8) || any(volbapar==9) || any(volbapar==14) || any(volbapar==15)
                    errordlg('Commission error of gravitational or disturbing tensor in the LNOF cannot be computed.',...
                        'Error in calculated parameters and output selection');
                    error('Commission error of gravitational or disturbing tensor in the LNOF cannot be computed.')
                end               
                
                %Import of the variance-covariance matrix file                           
                GGMcovname=get(findobj('tag','R'),'userdata');
                GGMcovadresar=get(findobj('tag','ell'),'userdata');

                %Error, if the input file has not been inputted
                if isempty(GGMcovname) 
                    errordlg('Please input a file with variance-covariance matrix of a geopotential model.',...
                        'Error in point type selection');
                    error('Please input a file with variance-covariance matrix of a geopotential model.')
                end

                %Update of the progress bar
                set(findobj('tag','hlasky'),'string',...
                        'Loading covariance matrix...','fontsize',8,...
                        'foregroundcolor','k'); drawnow;

                %Loading variance-covariance matrix
                if strcmp(GGMcovname(find(GGMcovname=='.'):end),'.mat') %If binary MAT file has been imported
                    covmat=importdata([GGMcovadresar,GGMcovname]);
                    nmaxGGMcov=max(covmat(:,2)); %Maximum degree of variance-covariance matrix
          
                    CSnm=covmat(:,1:3);
                    covmat=covmat(:,4:end);
                else %If ASCII file has been imported
                    fid=fopen([GGMcovadresar,GGMcovname]); 
                    covmat=fscanf(fid,'%f',[1 inf]);
                    fclose(fid);         

                    %Transformation of variance-covariance matrix from
                    %raw vector into the matrix
                    covmat=[0 0 0 0 0 0 covmat];

                    rows_GGM=length(covmat);
                    for i=0:2147483640
                        x=i^2-sum(1:(i-1));
                        if x==rows_GGM
                            rozmer=i;
                            break    	
                        end
                    end

                    korene=roots([1 2 -rozmer]);
                    nmaxGGMcov=korene(korene>0);

                    TM=triu(ones(nmaxGGMcov^2+2*nmaxGGMcov),0); %Lower triangular matrix with ones and zeros
                    TM(TM==1)=covmat;
                    covmat=TM'; %Lower triangular matrix with variances and covariances (raws 1:3 and columns 1:3 contain zeros)
                    
                    clear TM
                    
                    covmat=covmat(4:end,:);
  
                    if covmat(1,3)==0 && covmat(2,3)==1 %Sorted primarily according to degrees
                        errordlg('Wrong format of input variance-covariance matrix.',...
                            'Geopotential model and reference system selection');
                        error('Wrong format of input variance-covariance matrix.')
                    end
                    
                    CSnm=covmat(:,1:3);
                    covmat=covmat(:,4:end);
                end                
                                
                %Value of nmin and error messages
                nmin=str2double(get(findobj('tag','nmin'),'string'));
                if nmin<2 
                    errordlg('To compute commission error of a functional, the value of nmin must be at least 2.',...
                        'Error in geopotential model and reference system selection')
                    error('To compute commission error of a functional, the value of nmin must be at least 2.')
                elseif nmin>nmaxGGMcov;
                    errordlg('Value of nmin exceedes nmax value of the variance-covariance matrix of the GGM.',....
                        'Error in geopotential model and reference system selection')
                    error('Value of nmin exceedes nmax value of the variance-covariance matrix of the GGM.')
                end
                if isnan(nmin)==1
                        errordlg('Please input the nmin value.',...
                            'Error in geopotential model and reference system selection')
                        error('Please input the nmin value.')
                end

                nmin0=0; %Logical 0
                nmin1=0; %Logical 0
                
                %Value of nmax and error messages 
                if get(findobj('tag','use'),'value')==1
                    nmax=nmaxGGMcov;
                else
                    nmax=str2double(get(findobj('tag','nmax'),'string'));

                    if nmax>nmaxGGMcov
                        errordlg('Entered value of nmax exceedes nmax of the variance-covariance matrix.',...
                            'Error in geopotential model and reference system selection')
                        error('Entered value of nmax exceedes nmax of the variance-covariance matrix.')
                    elseif nmax<2
                        errordlg('Value of nmax must be at least 2.',...
                            'Error in geopotential model and reference system selection')
                        error('Value of nmax must be at least 2.')
                    elseif nmin>nmax
                        errordlg('Value of nmin cannot be larger than nmax value.',...
                            'Error in geopotential model and reference system selection')
                        error('Value of nmin cannot be larger than nmax value.')
                    end
                    if isnan(nmax)==1
                        errordlg('Please input the nmax value.',...
                            'Error in geopotential model and reference system selection')
                        error('Please input the nmax value.')
                    end
                end 
               
                %Deleting of variances and covariaces if nmax < nmaxGGM
                if nmax<nmaxGGMcov
                    zmazat=CSnm(:,2)>nmax;
                    covmat(zmazat,:)=[];
                    covmat(:,zmazat)=[];
                end

                clear CSnm zmazat

                [rows_cov cols_cov]=size(covmat);
                covmat=covmat'.*triu(ones(rows_cov),1)+covmat; %Creating full variances-covariance matrix from triangular matrix

                set(findobj('tag','hlasky'),'string',...
                        '','fontsize',8,'foregroundcolor','k'); drawnow;
  
                %Computation of commission error on a regular grid    
                if volbagridcheck==1      
                    
                    %Entered coordinates of the grid
                    fimin=str2double(get(findobj('tag','fimin'),'string')); 
                    fistep=str2double(get(findobj('tag','fistep'),'string'));
                    fimax=str2double(get(findobj('tag','fimax'),'string'));
                    lambdamin=str2double(get(findobj('tag','lambdamin'),'string')); 
                    lambdastep=str2double(get(findobj('tag','lambdastep'),'string'));
                    lambdamax=str2double(get(findobj('tag','lambdamax'),'string'));
                    h=str2double(get(findobj('tag','hgrid'),'string'));

                    %Check of the entered coordinates
                    if isnan(fimin) || isnan(fistep) || isnan(fimax) || isnan(lambdamin) || isnan(lambdastep) || isnan(lambdamax) || isnan(h)
                        errordlg('Entered grid is not correct.','Error in point type selection');
                        error('Entered grid is not correct.');       
                    end

                    if fimin>fimax
                        errordlg('Value of Lat. min must be smaller than the Lat. max value.',...
                            'Error in point type selection');
                        error('Value of Lat. min must be smaller than the Lat. max value.'); 
                    elseif fistep<=0
                        errordlg('Lat. step must be larger than zero.',...
                            'Error in point type selection');
                        error('Lat. step must be larger than zero.'); 
                    elseif lambdamin>lambdamax
                        errordlg('Value of Lon. min must be smaller than the Lon. max value.',...
                            'Error in point type selection');
                        error('Value of Lon. min must be smaller than the Lon. max value.');
                    elseif lambdastep<=0
                        errordlg('Lon. step must be larger than zero.',...
                            'Error in point type selection');
                        error('Lon. step must be larger than zero.'); 
                    end

                    if fimin>90 || fimin<-90
                        errordlg('Value of Lat. min must be within the interval <-90�,90�>.',...
                            'Error in point type selection');
                        error('Value of Lat. min must be within the interval <-90�,90�>.');
                    end
                    if fimax>90 || fimax<-90
                        errordlg('Value of Lat. max must be within the interval <-90�,90�>.',...
                            'Error in point type selection');
                        error('Value of Lat. max must be within the interval <-90�,90�>.');
                    end
                    if lambdamin>360 || lambdamin<-180
                        errordlg('Value of Lon. min must be within the interval <-180�,180�> or <0�,360�>.',...
                            'Error in point type selection');
                        error('Value of Lon. min must be within the interval <-180�,180�> or <0�,360�>.');
                    end
                    if lambdamax>360 || lambdamax<-180
                        errordlg('Value of Lon. max must be within the interval <-180�,180�> or <0�,360�>.',...
                            'Error in point type selection');
                        error('Value of Lon. max must be within the interval <-180�,180�> or <0�,360�>.');
                    end
                    if (lambdamax-lambdamin)>360
                        errordlg('Longitude must be in the range <-180�,180�> or <0�,360�>.',...
                            'Error in point type selection');
                        error('Longitude must be in the range <-180�,180�> or <0�,360�>.');
                    end
                    
                    %Vectors phi and lambda in one longitude and latitude
                    %parallel, respectively
                    fi=[fimin:fistep:fimax]';
                    lambda=[lambdamin:lambdastep:lambdamax]';  
                    [fi lambda]=meshgrid(fi,lambda);
                    [i_fi j_fi]=size(fi);
                    [i_lambda j_lambda]=size(lambda);
                    fi=fi'; 
                    fi=deg2rad(fi(:));
                    lambda=lambda'; 
                    lambda=deg2rad(lambda(:));  
                    
                    if coord==1
                        %Spherical radius
                        r=(R+h)*ones(length(fi),1);
                        hsph=h;
                        
                        %Spherical latitude
                        fiG=fi;  
                        
                        %Transform spherical latitude into the ellipsoidal
                        %latitude
                        [X Y Z]=sph2cart(0*fiG,fiG,r);
                        [fi , ~, h]=ecef2geodetic(X,Y,Z,[aEl eEl]');
                        clear X Y Z
                    elseif coord==0
                        %Transformation of the ellipsoidal coordinates into the
                        %cartesian coordinates
                        [X,Y,Z]=geodetic2ecef(fi,lambda,h,[aEl eEl]');
                        r=sqrt(X.*X+Y.*Y+Z.*Z); %Spherical radius
                        fiG=atan(Z./sqrt(X.*X+Y.*Y)); %Spherical latitude

                        clear X Y Z
                    end
                    
                    %If geoid/quasigeoid is to be computed
                    if any(volbapar==10) || any(volbapar==23)
                        if any(h~=0)
                            errordlg('In order to compute commission error of Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                                'Error in point type selection');
                            error('In order to compute commission error of Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                        end
                    end
                end
                
                %Computation of commission error in point-wise                            
                if volbadiskcheck==1               
                    %Ellipsoidal coordinates
                    fi=str2num(get(findobj('tag','fi'),'string'))';
                    lambda=str2num(get(findobj('tag','lambda'),'string'))';
                    h=str2num(get(findobj('tag','hdisk'),'string'))';
                end

                %Computation of commission error from imported data 
                if volbaloadcheck==1
                    %Import of data file containing ellipsoidal coordinates
                    loadname=get(findobj('tag','use'),'userdata');
                    loadadresar=get(findobj('tag','diskcheck'),'userdata');
                    Import=load([loadadresar,loadname],'-ascii');

                    %Ellipsoidal coordinates
                    fi=Import(:,1);
                    lambda=Import(:,2);
                    h=Import(:,3);                    
                    
                    volbadiskcheck=1;
                end

                if volbaloadcheck==1 || volbadiskcheck==1
                    
                    if any(fi>90) || any(fi<-90) 
                        errordlg('Values of Latitude must be within the interval <-90�,90�>.',...
                            'Error in point type selection');
                        error('Values of Latitude must be within the interval <-90�,90�>.');
                    end
                    if any(lambda>360) || any(lambda<-180)
                        errordlg('Values of Longitude min must be within the interval <-180�,180�> or <0�,360�>.',...
                            'Error in point type selection');
                        error('Values of Longitude min must be within the interval <-180�,180�> or <0�,360�>.');
                    end
                    
                    fi=deg2rad(fi(:));
                    lambda=deg2rad(lambda(:));  
                    
                    if coord==1 %Entered spherical coordinates                       
                        %Spherical radius
                        r=h;

                        %Spherical latitude
                        fiG=fi;  
                        
                        %Transform spherical latitude into the ellipsoidal
                        %latitude
                        [X Y Z]=sph2cart(lambda,fiG,r);
                        [fi , ~, h]=ecef2geodetic(X,Y,Z,[aEl eEl]');

                        clear X Y Z
                    elseif coord==0 %Entered ellipsoidal coordinates
                        %Trasformation of (fi, lambda, h) into (X, Y, Z)
                        [X,Y,Z]=geodetic2ecef(fi,lambda,h,[aEl eEl]');  
                        r=sqrt(X.*X+Y.*Y+Z.*Z); %Radius

                        %Spherical latitude
                        fiG=atan(Z./sqrt(X.*X+Y.*Y)); 
                        
                        clear X Y Z 
                    end
                end
 
                %=============================================================
                if volbadiskcheck==1                             

                    %Error message for inputted ellipsoidal coordinates
                    if isempty(fi) || isempty(lambda) || isempty(h)
                        errordlg('At least one ellipsoidal coordinate is empty or incorrectly entered.',...
                            'Error in point type selection');
                        error('At least one ellipsoidal coordinate is empty or incorrectly entered.');
                    end

                    if length(fi)~=length(lambda) || length(fi)~=length(h) || length(lambda)~=length(h)
                        errordlg('Ellipsoidal coordinates dimensions are not consistent or incorrectly entered.',...
                            'Error in point type selection')
                        error('Ellipsoidal coordinates dimensions are not consistent or incorrectly entered.')     
                    end               
                    
                    %If geoid/quasigeoid is to be computed
                    if any(volbapar==10) || any(volbapar==23)
                        if any(h~=0)
                            errordlg('In order to compute commission error of Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                                'Error in point type selection');
                            error('In order to compute commission error of Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                        end
                    end
                end
                                
                %Initialization of the matrices and vectors for the computation of fnALFs
                Pnm=zeros(length(fi),nmax+1);
                q=(R./r);
                q2=(R./r).^2;
                u=cos(fiG);
                t=sin(fiG);
                
                %Initialization of the matrices and vectors for the 
                %computation of the first-order derivatives of fnALFs
                if any(volbapar==3) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16)  || any(volbapar==20)
                    dALFs=1;
                    dPnm=zeros(length(fi),nmax+1);
                    qu=q./u;
                    tu=t./u;
                    
                    %Treatment of the dPnm singularity
                    singdPnm=fi==pi/2 | fi==-pi/2;
                else
                    dALFs=0;
                end 
                
                %Initialization of the matrices and vectors for the 
                %computation of the second-order derivatives of fnALFs
                if any(volbapar==6) || any(volbapar==12)
                    ddALFs=1;
                    ddPnm=zeros(length(fi),nmax+1);
                    
                    %Treatment of the ddPnm singularity
                    singddPnm=fi==pi/2 | fi==-pi/2;
                else
                    ddALFs=0;
                end  
                
                z=((nmax+1)*(nmax+2)-6)/2-(nmax-1);
                k=z;
                
                %Status line
                progressbar=findobj('tag','hlasky');
                
                %% Summation over m
                for m=nmax:-1:0
                    
                    %Update of the progress bar
                    if rem(m,10)==0
                        set(progressbar,'string',...
                            sprintf('Progress: m = %5.0d',m),...
                            'fontsize',8); drawnow;
                    end
                    
                    if m==1
                        z=z+1;
                    end
                    
                    %Standard forward column method
                    if m==0
                        Pnm(:,1)=1;
                    elseif m==1                    
                        Pnm(:,1)=sqrt(3)*u.*q;  
                    elseif m>1                            
                        i=2*(2:m);
                        i1=sqrt((i+ones(size(i)))./i);
                        Pnm(:,1)=u.^m*sqrt(3)*prod(i1).*q.^m;
                    end

                    if m==nmax
                    elseif m<=(nmax-1)
                        n=m+1;
                        anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                        Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                    end

                    if m<(nmax-1)
                        j=3;
                        for n=m+2:nmax                            
                            anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                            bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                            Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                            j=j+1;
                        end
                    end               
                      
                    %% Computation of the first-order derivatives of modified fnALFs
                    if dALFs==1  
                        if volbaALFs==1 || volbaALFs==3
                            enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                            if m==0 %Zonal modified dALFs
                                dPnm(:,1)=0.*u;
                                dPnm(:,2)=sqrt(3)*u.*q;
                                dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                            elseif m==nmax %Sectorial modified dALFs
                                dPnm(:,1)=-m*(tu).*Pnm(:,1);
                            else
                                dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne ALFs
                            end

                        elseif volbaALFs==2
                            enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                            if m==0 %Zonal modified dALFs
                                dPnm(:,1)=0.*u*1e-280;
                                dPnm(:,2)=sqrt(3)*u.*q*1e-280;
                                dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                            elseif m==nmax
                                dPnm(:,1)=-m*(tu).*Pnm(:,1); %Sectorial modified dALFs
                            else
                                dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne dALFs
                            end
                        end                 

                        %Treatment of the dPnm singularity
                        dPnm(singdPnm,:)=0;
                        
                        
                        %Computation of the second-order derivatives of modified fnALFs
                        if ddALFs==1
                            if m==0 %Zonal modified ddALFs
                                ddPnm=bsxfun(@times,tu,dPnm)-bsxfun(@times,(0:nmax).*((0:nmax)+1),Pnm);
                            else                                  
                                ddPnm(:,1:end-m)=bsxfun(@times,tu,dPnm(:,1:end-m))+bsxfun(@times,m^2./u.^2,Pnm(:,1:end-m))-bsxfun(@times,(m:nmax).*((m:nmax)+1),Pnm(:,1:end-m));
                            end
                                                                
                            %Treatment of the ddPnm singularity
                            ddPnm(singddPnm,:)=0;                    
                        end
                        
                    end
                    
                    %If nmin>2
                    if nmin~=2
                        Pnm(:,1:(nmin-m))=0;
                        
                        if dALFs==1
                            dPnm(:,1:(nmin-m))=0;
                            if ddALFs==1
                                ddPnm(:,1:(nmin-m))=0;
                            end
                        end
                    end
             
                    cosla=cos(m*lambda);
                    sinla=sin(m*lambda);

                    %% Loop for 1:NF (number of computing functionals)
                    for i=1:pocetpar                   
                        if volbapar(i)==1   
                        elseif volbapar(i)==2 %Deflection of the vertical eta                           

                            if m==nmax
                                ACeta=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASeta=ACeta;
                                Aeta=zeros(length(fiG),cols_cov-(nmax-1));
                            end

                            if m==0
                                Aeta(:,1:2:end)=ACeta;
                                Aeta(:,2:2:end)=ASeta;

                                Aeta=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla) Aeta];
                            elseif m==1
                                ACeta(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla);
                                ASeta(:,z:k)=-bsxfun(@times,m*Pnm(:,(3-m):(end-m)),cosla); 
                            else
                                ACeta(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),sinla);
                                ASeta(:,z:k)=-bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),cosla);
                            end

                            if m==0
                                clear ACeta ASeta
                            end

                        elseif volbapar(i)==3 %Deflection of the vertical xi

                            if m==nmax
                                ACksi=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASksi=ACksi;
                                Aksi=zeros(length(fiG),cols_cov-(nmax-1));
                            end

                            if m==0
                                Aksi(:,1:2:end)=ACksi;
                                Aksi(:,2:2:end)=ASksi;

                                Aksi=[bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla) Aksi];
                            elseif m==1
                                ACksi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla);
                                ASksi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACksi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),cosla);
                                ASksi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACksi ASksi
                            end

                        elseif volbapar(i)==4 %Deflection of the vertical Theta

                            if m==nmax
                                ACTeta=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASTeta=ACTeta;
                                ACTksi=ACTeta;
                                ASTksi=ACTeta;
                                ATeta=zeros(length(fiG),cols_cov-(nmax-1));
                                ATksi=ATeta;
                            end

                            if m==0
                                ATeta(:,1:2:end)=ACTeta;
                                ATeta(:,2:2:end)=ASTeta;
                                ATksi(:,1:2:end)=ACTksi;
                                ATksi(:,2:2:end)=ASTksi;

                                ATksi=[bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla) ATksi];
                                ATeta=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla) ATeta];
                            elseif m==1
                                ACTksi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla);
                                ASTksi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),sinla); 
                                
                                ACTeta(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla);
                                ASTeta(:,z:k)=-bsxfun(@times,m*Pnm(:,(3-m):(end-m)),cosla); 
                            else
                                ACTksi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),cosla);
                                ASTksi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACTeta(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),sinla);
                                ASTeta(:,z:k)=-bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),cosla);
                            end

                            if m==0
                                clear ACTeta ASTeta ACTksi ASTksi
                            end

                        elseif volbapar(i)==5 %Disturbing potential

                            if m==nmax
                                ACT=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                AST=ACT;
                                AT=zeros(length(fiG),cols_cov-(nmax-1));
                            end

                            if m==0
                                AT(:,1:2:end)=ACT;
                                AT(:,2:2:end)=AST;

                                AT=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AT];
                            elseif m==1
                                ACT(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                AST(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACT(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                AST(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACT AST
                            end
                            
                        elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                            
                            if m==nmax
                                ACTrr=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASTrr=ACTrr;
                                ACTfifi=ACTrr;
                                ASTfifi=ACTrr;
                                ATrr=zeros(length(fiG),cols_cov-(nmax-1));
                                ATfifi=ATrr;
                                ACTll=ACTrr;
                                ASTll=ACTrr;
                                ATll=ATrr;
                                ampl_Trr=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Trr(:,n-1)=(n+1)*(n+2);
                                end
                            end

                            if m==0
                                ATrr(:,1:2:end)=ACTrr;
                                ATrr(:,2:2:end)=ASTrr;

                                ATrr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Trr,cosla) ATrr];
                                
                                ATfifi(:,1:2:end)=ACTfifi;
                                ATfifi(:,2:2:end)=ASTfifi;
                                
                                ATfifi=[bsxfun(@times,ddPnm(:,(3-m):(end-m)),cosla) ATfifi];
                                
                                ATll(:,1:2:end)=ACTll;
                                ATll(:,2:2:end)=ASTll;

                                ATll=[bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),cosla) ATll];
                            elseif m==1
                                ACTrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Trr,cosla);
                                ASTrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Trr,sinla); 
                                
                                ACTfifi(:,z:k)=bsxfun(@times,ddPnm(:,(3-m):(end-m)),cosla);
                                ASTfifi(:,z:k)=bsxfun(@times,ddPnm(:,(3-m):(end-m)),sinla); 
                                
                                ACTll(:,z:k)=bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),cosla);
                                ASTll(:,z:k)=bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACTrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Trr(:,(m-1):end),cosla);
                                ASTrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Trr(:,(m-1):end),sinla);
                                
                                ACTfifi(:,z:k)=bsxfun(@times,ddPnm(:,1:(nmax+1-m)),cosla);
                                ASTfifi(:,z:k)=bsxfun(@times,ddPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACTll(:,z:k)=bsxfun(@times,m^2*Pnm(:,1:(nmax+1-m)),cosla);
                                ASTll(:,z:k)=bsxfun(@times,m^2*Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ampl_Trr ACTrr ASTrr ACTfifi ...
                                    ASTfifi ACTll ASTll
                            end
                            
                        elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                            
                            if m==nmax
                                ACTrfi=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASTrfi=ACTrfi;
                                ATrfi=zeros(length(fiG),cols_cov-(nmax-1));
                                ACTrl=ACTrfi;
                                ASTrl=ACTrfi;
                                ATrl=ATrfi;
                                ACTfil=ACTrfi;
                                ASTfil=ACTrfi;
                                ATfil=ATrfi;
                                ampl_Tr=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Tr(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                ATrfi(:,1:2:end)=ACTrfi;
                                ATrfi(:,2:2:end)=ASTrfi;

                                ATrfi=[bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Tr,cosla) ATrfi];
                                
                                ATrl(:,1:2:end)=ACTrl;
                                ATrl(:,2:2:end)=ASTrl;

                                ATrl=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Tr,sinla) ATrl];
                                
                                ATfil(:,1:2:end)=ACTfil;
                                ATfil(:,2:2:end)=ASTfil;

                                ATfil=[bsxfun(@times,m*dPnm(:,(3-m):(end-m)),sinla) ATfil];
                            elseif m==1
                                ACTrfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Tr,cosla);
                                ASTrfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Tr,sinla); 

                                ACTrl(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Tr,sinla);
                                ASTrl(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Tr,cosla); 
                                
                                ACTfil(:,z:k)=bsxfun(@times,m*dPnm(:,(3-m):(end-m)),sinla);
                                ASTfil(:,z:k)=bsxfun(@times,m*dPnm(:,(3-m):(end-m)),cosla); 
                            else
                                ACTrfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)).*ampl_Tr(:,(m-1):end),cosla);
                                ASTrfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)).*ampl_Tr(:,(m-1):end),sinla);

                                ACTrl(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)).*ampl_Tr(:,(m-1):end),sinla);
                                ASTrl(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)).*ampl_Tr(:,(m-1):end),cosla);
                                
                                ACTfil(:,z:k)=bsxfun(@times,m*dPnm(:,1:(nmax+1-m)),sinla);
                                ASTfil(:,z:k)=bsxfun(@times,m*dPnm(:,1:(nmax+1-m)),cosla);
                            end

                            if m==0
                                clear ampl_Trfi ACTrfi ASTrfi ACTrl ASTrl ... 
                                    ACTfil ASTfil ampl_Tr
                            end
                            
                        elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                        elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                        elseif volbapar(i)==10 %Geoid undulation

                            if m==nmax
                                ACN=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASN=ACN;
                                AN=zeros(length(fiG),cols_cov-(nmax-1));
                            end

                            if m==0
                                AN(:,1:2:end)=ACN;
                                AN(:,2:2:end)=ASN;

                                AN=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AN];
                            elseif m==1
                                ACN(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASN(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACN(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASN(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACN ASN
                            end
                            
                        elseif volbapar(i)==11 %Gravitational potential

                            if m==nmax
                                ACV=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASV=ACV;
                                AV=zeros(length(fiG),cols_cov-(nmax-1));
                            end

                            if m==0
                                AV(:,1:2:end)=ACV;
                                AV(:,2:2:end)=ASV;

                                AV=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AV];
                            elseif m==1
                                ACV(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASV(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACV(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASV(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACV ASV
                            end

                        elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                            
                            if m==nmax
                                ACVrr=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASVrr=ACVrr;
                                ACVfifi=ACVrr;
                                ASVfifi=ACVrr;
                                AVrr=zeros(length(fiG),cols_cov-(nmax-1));
                                AVfifi=AVrr;
                                ACVll=ACVrr;
                                ASVll=ACVrr;
                                AVll=AVrr;
                                ampl_Vrr=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Vrr(:,n-1)=(n+1)*(n+2);
                                end
                            end

                            if m==0
                                AVrr(:,1:2:end)=ACVrr;
                                AVrr(:,2:2:end)=ASVrr;

                                AVrr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Vrr,cosla) AVrr];
                                
                                AVfifi(:,1:2:end)=ACVfifi;
                                AVfifi(:,2:2:end)=ASVfifi;
                                
                                AVfifi=[bsxfun(@times,ddPnm(:,(3-m):(end-m)),cosla) AVfifi];
                                
                                AVll(:,1:2:end)=ACVll;
                                AVll(:,2:2:end)=ASVll;

                                AVll=[bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),cosla) AVll];
                            elseif m==1
                                ACVrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Vrr,cosla);
                                ASVrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Vrr,sinla); 
                                
                                ACVfifi(:,z:k)=bsxfun(@times,ddPnm(:,(3-m):(end-m)),cosla);
                                ASVfifi(:,z:k)=bsxfun(@times,ddPnm(:,(3-m):(end-m)),sinla); 
                                
                                ACVll(:,z:k)=bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),cosla);
                                ASVll(:,z:k)=bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACVrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Vrr(:,(m-1):end),cosla);
                                ASVrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Vrr(:,(m-1):end),sinla);
                                
                                ACVfifi(:,z:k)=bsxfun(@times,ddPnm(:,1:(nmax+1-m)),cosla);
                                ASVfifi(:,z:k)=bsxfun(@times,ddPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACVll(:,z:k)=bsxfun(@times,m^2*Pnm(:,1:(nmax+1-m)),cosla);
                                ASVll(:,z:k)=bsxfun(@times,m^2*Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ampl_Vrr ACVrr ASVrr ACVfifi ...
                                    ASVfifi ACVll ASVll
                            end
                            
                        elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                            
                            if m==nmax
                                ACVrfi=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASVrfi=ACVrfi;
                                AVrfi=zeros(length(fiG),cols_cov-(nmax-1));
                                ACVrl=ACVrfi;
                                ASVrl=ACVrfi;
                                AVrl=AVrfi;
                                ACVfil=ACVrfi;
                                ASVfil=ACVrfi;
                                AVfil=AVrfi;
                                ampl_Vr=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Vr(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                AVrfi(:,1:2:end)=ACVrfi;
                                AVrfi(:,2:2:end)=ASVrfi;

                                AVrfi=[bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Vr,cosla) AVrfi];
                                
                                AVrl(:,1:2:end)=ACVrl;
                                AVrl(:,2:2:end)=ASVrl;

                                AVrl=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Vr,sinla) AVrl];
                                
                                AVfil(:,1:2:end)=ACVfil;
                                AVfil(:,2:2:end)=ASVfil;

                                AVfil=[bsxfun(@times,m*dPnm(:,(3-m):(end-m)),sinla) AVfil];
                            elseif m==1
                                ACVrfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Vr,cosla);
                                ASVrfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Vr,sinla); 

                                ACVrl(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Vr,sinla);
                                ASVrl(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Vr,cosla); 
                                
                                ACVfil(:,z:k)=bsxfun(@times,m*dPnm(:,(3-m):(end-m)),sinla);
                                ASVfil(:,z:k)=bsxfun(@times,m*dPnm(:,(3-m):(end-m)),cosla); 
                            else
                                ACVrfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)).*ampl_Vr(:,(m-1):end),cosla);
                                ASVrfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)).*ampl_Vr(:,(m-1):end),sinla);

                                ACVrl(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)).*ampl_Vr(:,(m-1):end),sinla);
                                ASVrl(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)).*ampl_Vr(:,(m-1):end),cosla);
                                
                                ACVfil(:,z:k)=bsxfun(@times,m*dPnm(:,1:(nmax+1-m)),sinla);
                                ASVfil(:,z:k)=bsxfun(@times,m*dPnm(:,1:(nmax+1-m)),cosla);
                            end

                            if m==0
                                clear ampl_Vrfi ACVrfi ASVrfi ACVrl ...
                                    ASVrl ACVfil ASVfil ampl_Vr
                            end
                            
                        elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                        elseif volbapar(i)==15 %Vxy_Vxz_Vyz
                        elseif volbapar(i)==16 %Gravity  

                            if m==nmax
                                ACWr=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASWr=ACWr;
                                ACWfi=ACWr;
                                ASWfi=ACWr;
                                ACWlambda=ACWr;
                                ASWlambda=ACWr;
                                AWr=zeros(length(fiG),cols_cov-(nmax-1));
                                AWfi=AWr;
                                AWlambda=AWr;
                                ampl_Wr=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Wr(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                AWr(:,1:2:end)=ACWr;
                                AWr(:,2:2:end)=ASWr;
                                AWfi(:,1:2:end)=ACWfi;
                                AWfi(:,2:2:end)=ASWfi;
                                AWlambda(:,1:2:end)=ACWlambda;
                                AWlambda(:,2:2:end)=ASWlambda;

                                AWlambda=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla) AWlambda];
                                AWfi=[bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla) AWfi];
                                AWr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr,cosla) AWr];
                            elseif m==1
                                ACWlambda(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla);
                                ASWlambda(:,z:k)=-bsxfun(@times,m*Pnm(:,(3-m):(end-m)),cosla); 
                                
                                ACWfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla);
                                ASWfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),sinla);
                                
                                ACWr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr,cosla);
                                ASWr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr,sinla); 
                            else
                                ACWlambda(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),sinla);
                                ASWlambda(:,z:k)=-bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),cosla);
                                
                                ACWfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),cosla);
                                ASWfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACWr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wr(:,(m-1):end),cosla);
                                ASWr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wr(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_Wr ACWr ASWr ACWfi ASWfi ...
                                    ACWlambda ASWlambda
                            end
                                                        
                        elseif volbapar(i)==17 %Gravity sa

                            if m==nmax
                                ACg_sa=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASg_sa=ACg_sa;
                                Ag_sa=zeros(length(fiG),cols_cov-(nmax-1));
                                ampl_g_sa=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_g_sa(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                Ag_sa(:,1:2:end)=ACg_sa;
                                Ag_sa(:,2:2:end)=ASg_sa;

                                Ag_sa=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_g_sa,cosla) Ag_sa];
                            elseif m==1
                                ACg_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_g_sa,cosla);
                                ASg_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_g_sa,sinla); 
                            else
                                ACg_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_g_sa(:,(m-1):end),cosla);
                                ASg_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_g_sa(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_g_sa ACg_sa ASg_sa
                            end

                        elseif volbapar(i)==18 %Gravity potential                            

                            if m==nmax
                                ACW=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASW=ACW;
                                AW=zeros(length(fiG),cols_cov-(nmax-1));
                            end

                            if m==0
                                AW(:,1:2:end)=ACW;
                                AW(:,2:2:end)=ASW;

                                AW=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AW];
                            elseif m==1
                                ACW(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASW(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACW(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASW(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACW ASW
                            end

                        elseif volbapar(i)==19 %Gravity anomaly sa

                            if m==nmax
                                ACanomalia_sa=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASanomalia_sa=ACanomalia_sa;
                                Aanomalia_sa=zeros(length(fiG),cols_cov-(nmax-1));
                                ampl_anomalia_sa=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_anomalia_sa(:,n-1)=(n-1);
                                end
                            end

                            if m==0
                                Aanomalia_sa(:,1:2:end)=ACanomalia_sa;
                                Aanomalia_sa(:,2:2:end)=ASanomalia_sa;

                                Aanomalia_sa=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_anomalia_sa,cosla) Aanomalia_sa];
                            elseif m==1
                                ACanomalia_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_anomalia_sa,cosla);
                                ASanomalia_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_anomalia_sa,sinla); 
                            else
                                ACanomalia_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_anomalia_sa(:,(m-1):end),cosla);
                                ASanomalia_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_anomalia_sa(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_anomalia_sa ACanomalia_sa ...
                                    ASanomalia_sa
                            end

                        elseif volbapar(i)==20 %Gravity disturbance

                            if m==nmax
                                ACWr_porucha=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASWr_porucha=ACWr_porucha;
                                ACWfi_porucha=ACWr_porucha;
                                ASWfi_porucha=ACWr_porucha;
                                ACWlambda_porucha=ACWr_porucha;
                                ASWlambda_porucha=ACWr_porucha;
                                AWr_porucha=zeros(length(fiG),cols_cov-(nmax-1));
                                AWfi_porucha=AWr_porucha;
                                AWlambda_porucha=AWr_porucha;
                                ampl_Wr_porucha=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Wr_porucha(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                AWr_porucha(:,1:2:end)=ACWr_porucha;
                                AWr_porucha(:,2:2:end)=ASWr_porucha;
                                AWfi_porucha(:,1:2:end)=ACWfi_porucha;
                                AWfi_porucha(:,2:2:end)=ASWfi_porucha;
                                AWlambda_porucha(:,1:2:end)=ACWlambda_porucha;
                                AWlambda_porucha(:,2:2:end)=ASWlambda_porucha;

                                AWlambda_porucha=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla) AWlambda_porucha];
                                AWfi_porucha=[bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla) AWfi_porucha];
                                AWr_porucha=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr_porucha,cosla) AWr_porucha];
                            elseif m==1
                                ACWlambda_porucha(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla);
                                ASWlambda_porucha(:,z:k)=-bsxfun(@times,m*Pnm(:,(3-m):(end-m)),cosla); 
                                
                                ACWfi_porucha(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla);
                                ASWfi_porucha(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),sinla);
                                
                                ACWr_porucha(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr_porucha,cosla);
                                ASWr_porucha(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr_porucha,sinla); 
                            else
                                ACWlambda_porucha(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),sinla);
                                ASWlambda_porucha(:,z:k)=-bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),cosla);
                                
                                ACWfi_porucha(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),cosla);
                                ASWfi_porucha(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACWr_porucha(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wr_porucha(:,(m-1):end),cosla);
                                ASWr_porucha(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wr_porucha(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_Wr_porucha ACWr_porucha ...
                                    ASWr_porucha ACWfi_porucha ...
                                    ASWfi_porucha ACWlambda_porucha ASWlambda
                            end

                        elseif volbapar(i)==21 %Gravity disturbance sa

                            if m==nmax
                                ACporucha_sa=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASporucha_sa=ACporucha_sa;
                                Aporucha_sa=zeros(length(fiG),cols_cov-(nmax-1));
                                ampl_porucha_sa=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_porucha_sa(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                Aporucha_sa(:,1:2:end)=ACporucha_sa;
                                Aporucha_sa(:,2:2:end)=ASporucha_sa;

                                Aporucha_sa=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_porucha_sa,cosla) Aporucha_sa];
                            elseif m==1
                                ACporucha_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_porucha_sa,cosla);
                                ASporucha_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_porucha_sa,sinla); 
                            else
                                ACporucha_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_porucha_sa(:,(m-1):end),cosla);
                                ASporucha_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_porucha_sa(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_porucha_sa ACporucha_sa ASporucha_sa
                            end

                        elseif volbapar(i)==22 %Height anomaly Ell

                            if m==nmax
                                ACzetaEl=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASzetaEl=ACzetaEl;
                                AzetaEl=zeros(length(fiG),cols_cov-(nmax-1));
                            end

                            if m==0
                                AzetaEl(:,1:2:end)=ACzetaEl;
                                AzetaEl(:,2:2:end)=ASzetaEl;

                                AzetaEl=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AzetaEl];
                            elseif m==1
                                ACzetaEl(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASzetaEl(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACzetaEl(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASzetaEl(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACzetaEl ASzetaEl
                            end

                        elseif volbapar(i)==23 %Height anomaly
                            
                            if m==nmax
                                ACzeta=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASzeta=ACzeta;
                                Azeta=zeros(length(fiG),cols_cov-(nmax-1));
                            end

                            if m==0
                                Azeta(:,1:2:end)=ACzeta;
                                Azeta(:,2:2:end)=ASzeta;

                                Azeta=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) Azeta];
                            elseif m==1
                                ACzeta(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASzeta(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACzeta(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASzeta(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACzeta ASzeta
                            end
                            
                        elseif volbapar(i)==24 %Second radial derivative of disturbing potential

                            if m==nmax
                                ACT_rr=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                AST_rr=ACT_rr;
                                AT_rr=zeros(length(fiG),cols_cov-(nmax-1));
                                ampl_T_rr=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_T_rr(:,n-1)=(n+1)*(n+2);
                                end
                            end

                            if m==0
                                AT_rr(:,1:2:end)=ACT_rr;
                                AT_rr(:,2:2:end)=AST_rr;

                                AT_rr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_T_rr,cosla) AT_rr];
                            elseif m==1
                                ACT_rr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_T_rr,cosla);
                                AST_rr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_T_rr,sinla); 
                            else
                                ACT_rr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_T_rr(:,(m-1):end),cosla);
                                AST_rr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_T_rr(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_T_rr ACT_rr AST_rr
                            end

                        elseif volbapar(i)==25 %Second radial derivative of gravity potential

                            if m==nmax
                                ACWrr=zeros(length(fiG),((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASWrr=ACWrr;
                                AWrr=zeros(length(fiG),cols_cov-(nmax-1));
                                ampl_Wrr=zeros(length(fiG),nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Wrr(:,n-1)=(n+1)*(n+2);
                                end
                            end

                            if m==0
                                AWrr(:,1:2:end)=ACWrr;
                                AWrr(:,2:2:end)=ASWrr;

                                AWrr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wrr,cosla) AWrr];
                            elseif m==1
                                ACWrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wrr,cosla);
                                ASWrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wrr,sinla); 
                            else
                                ACWrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wrr(:,(m-1):end),cosla);
                                ASWrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wrr(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_Wrr ACWrr ASWrr
                            end
                        end
                    end
                    
                    if m==0
                    elseif m==1
                        z=z-nmax+m-1;
                        k=z+nmax-m; 
                    else
                        z=z-nmax+m-2;
                        k=z+nmax-m+1; 
                    end
                end
                
                %Update of progress bar
                set(progressbar,'string','Progress: Matrix multiplications...','fontsize',8); drawnow;
                
                clear Pnm dPnm ddPnm t u q q2 tu qu sinla cosla enm ...
                    singdPnm singddPnm
                
                %Computation of the normal gravity for eta, xi, Theta,
                %geoid undulation, zeta El, zeta
                if any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==10) || any(volbapar==22) || any(volbapar==23)                        
                    bEl=aEl*sqrt(1-eEl^2);
                    EEl=sqrt(aEl^2-bEl^2);
                    
                    if volbagridcheck==1
                        %Computation of the ellipsoidal harmonic coordinates
                        [X,Y,Z]=geodetic2ecef(fi,0*zeros(length(fi),1),h,[aEl eEl]');
                        ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                        betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                        clear X Y Z
                    elseif volbadiskcheck==1
                        %Computation of the ellipsoidal harmonic coordinates
                        [X,Y,Z]=geodetic2ecef(fi,lambda,h,[aEl eEl]');
                        ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                        betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                        clear X Y Z
                    end

                    wgama=sqrt((ugama.^2+EEl^2*sin(betagama).^2)./(ugama.^2+EEl^2));
                    qgama=1/2*((1+(3*ugama.^2)./EEl^2).*atan(EEl./ugama)-3*ugama./EEl);
                    qgama_=3*(1+(ugama.^2)./EEl^2).*(1-ugama./EEl.*atan(EEl./ugama))-1;
                    qgama0=1/2*((1+(3*bEl^2)/EEl^2)*atan(EEl/bEl)-3*bEl/EEl);
                        
                    gamau=-1./wgama.*(GMEl./(ugama.^2+EEl.^2)+(omegaEl.^2.*aEl.^2.*EEl)./(ugama.^2+EEl.^2).*(qgama_./qgama0).*(1./2.*sin(betagama).^2-1./6)-omegaEl.^2.*ugama.*cos(betagama).^2);
                    gamabeta=-1./wgama.*(-(omegaEl.^2.*aEl.^2.*qgama)./(sqrt(ugama.^2+EEl.^2).*qgama0)+omegaEl.^2.*sqrt(ugama.^2+EEl.^2)).*sin(betagama).*cos(betagama);
                        
                    clear ugama betagama wgama qgama qgama_ qgama0
                       
                    gamaP=sqrt(gamau.^2+gamabeta.^2);
                        
                    clear gamau gamabeta
                end
                
                %% Final computation of commission errors of the functionals             
                for i=1:pocetpar
                    if volbapar(i)==1                
                    elseif volbapar(i)==2 %Deflection of the vertical eta

                        Pg=sum((Aeta*covmat).*Aeta,2);
                        clear Aeta
                        Pg=(GM./(r.^2.*gamaP.*cos(fiG)).*sqrt(Pg))*(180/pi)*3600;
                        Pg(fi==pi/2 | fi==-pi/2)=0;
                        
                    elseif volbapar(i)==3 %Deflection of the vertical xi

                        Pg=sum((Aksi*covmat).*Aksi,2);
                        clear Aksi
                        Pg=(GM./(r.^2.*gamaP).*sqrt(Pg))*(180/pi)*3600;
                        
                    elseif volbapar(i)==4 %Deflection of the vertical Theta

                        eta=sum((ATeta*covmat).*ATeta,2);
                        clear ATeta
                       
                        ksi=sum((ATksi*covmat).*ATksi,2);
                        clear ATksi
                        
                        eta=(GM./(r.^2.*gamaP.*cos(fiG)).*sqrt(eta))*(180/pi)*3600;
                        eta(fi==pi/2 | fi==-pi/2)=0;
                        ksi=(GM./(r.^2.*gamaP).*sqrt(ksi))*(180/pi)*3600;

                        Pg=sqrt(eta.^2+ksi.^2);
                       
                        clear eta ksi
                    elseif volbapar(i)==5 %Disturbing potential                        
     
                        Pg=sum((AT*covmat).*AT,2);
                        clear AT
                        Pg=(GM./r.*sqrt(Pg));
                        
                    elseif volbapar(i)==6 %Disturbing tensor Trr_Tp_Tll
                        
                        Trr=sum((ATrr*covmat).*ATrr,2);
                        clear ATrr
                        Trr=(GM./r.^3.*sqrt(Trr))*10^9;
                        
                        Tfifi=sum((ATfifi*covmat).*ATfifi,2);
                        clear ATfifi
                        Tfifi=(GM./r.^3.*sqrt(Tfifi))*10^9;
                        
                        Tll=sum((ATll*covmat).*ATll,2);
                        clear ATll
                        Tll=(GM./(r.^3.*cos(fiG).^2).*sqrt(Tll))*10^9;
                        Tll(fi>deg2rad(89.9) | fi<deg2rad(-89.9))=0;
                        
                        Pg=Trr;
                        clear Trr
                        Pg=[Pg Tfifi];
                        clear Tfifi
                        Pg=[Pg Tll];
                        clear Tll
                        
                    elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                        
                        Trfi=sum((ATrfi*covmat).*ATrfi,2);
                        clear ATrfi
                        Trfi=(GM./r.^3.*sqrt(Trfi))*10^9;
                        
                        Trl=sum((ATrl*covmat).*ATrl,2);
                        clear ATrl
                        Trl=(GM./(r.^3.*cos(fiG)).*sqrt(Trl))*10^9;
                        
                        Tfil=sum((ATfil*covmat).*ATfil,2);
                        clear ATfil
                        Tfil=(GM./(r.^3.*cos(fiG)).*sqrt(Tfil))*10^9;
                        Tfil(fi>deg2rad(89.5) | fi<deg2rad(-89.5),:)=0;
                        
                        Pg=Trfi;
                        clear Trfi
                        Pg=[Pg Trl];
                        clear Trl
                        Pg=[Pg Tfil];
                        clear Tfil
                        
                    elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                    elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                    elseif volbapar(i)==10 %Geoid undulation
                        
                        Pg=sum((AN*covmat).*AN,2);
                        clear AN
                        Pg=(GM./(r.*gamaP).*sqrt(Pg));
                        
                    elseif volbapar(i)==11 %Gravitational potential
                        
                        Pg=sum((AV*covmat).*AV,2);
                        clear AV
                        Pg=(GM./r.*sqrt(Pg));
                        
                    elseif volbapar(i)==12 %Gravitational tensor Vrr_Vp_Vll
                        
                        Vrr=sum((AVrr*covmat).*AVrr,2);
                        clear AVrr
                        Vrr=(GM./r.^3.*sqrt(Vrr))*10^9;
                        
                        Vfifi=sum((AVfifi*covmat).*AVfifi,2);
                        clear AVfifi
                        Vfifi=(GM./r.^3.*sqrt(Vfifi))*10^9;
                        
                        Vll=sum((AVll*covmat).*AVll,2);
                        clear AVll
                        Vll=(GM./(r.^3.*cos(fiG).^2).*sqrt(Vll))*10^9;
                        Vll(fi>deg2rad(89.9) | fi<deg2rad(-89.9))=0;
                        
                        Pg=Vrr;
                        clear Vrr
                        Pg=[Pg Vfifi];
                        clear Vfifi
                        Pg=[Pg Vll];
                        clear Vll
                        
                    elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                        
                        Vrfi=sum((AVrfi*covmat).*AVrfi,2);
                        clear AVrfi
                        Vrfi=(GM./r.^3.*sqrt(Vrfi))*10^9;
                        
                        Vrl=sum((AVrl*covmat).*AVrl,2);
                        clear AVrl
                        Vrl=(GM./(r.^3.*cos(fiG)).*sqrt(Vrl))*10^9;
                        
                        Vfil=sum((AVfil*covmat).*AVfil,2);
                        clear AVfil
                        Vfil=(GM./(r.^3.*cos(fiG)).*sqrt(Vfil))*10^9;
                        Vfil(fi>deg2rad(89.5) | fi<deg2rad(-89.5),:)=0;
                        
                        Pg=Vrfi;
                        clear Vrfi
                        Pg=[Pg Vrl];
                        clear Vrl
                        Pg=[Pg Vfil];
                        clear Vfil
                        
                    elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                    elseif volbapar(i)==15 %Gravitational tensor Vxy_Vxz_Vyz
                    elseif volbapar(i)==16 %Gravity
                        
                        Wr=sum((AWr*covmat).*AWr,2);
                        clear AWr
                        Wr=(GM./r.^2.*sqrt(Wr));
                        
                        Wfi=sum((AWfi*covmat).*AWfi,2);
                        clear AWfi
                        Wfi=(GM./r.^2.*sqrt(Wfi));
                        
                        Wlambda=sum((AWlambda*covmat).*AWlambda,2);
                        clear AWlambda
                        Wlambda=(GM./(r.^2.*cos(fiG)).*sqrt(Wlambda));

                        Pg=sqrt(Wr.^2+Wlambda.^2+Wfi.^2)*10^5;

                        clear Wr Wfi Wlambda
                    elseif volbapar(i)==17 %Gravity sa
                        
                        Pg=sum((Ag_sa*covmat).*Ag_sa,2);
                        clear Ag_sa
                        Pg=(GM./r.^2.*sqrt(Pg))*10^5;
                        
                    elseif volbapar(i)==18 %Gravity potential
                        
                        Pg=sum((AW*covmat).*AW,2);
                        clear AW
                        Pg=(GM./r.*sqrt(Pg));
                        
                    elseif volbapar(i)==19 %Gravity anomaly sa
                        
                        Pg=sum((Aanomalia_sa*covmat).*Aanomalia_sa,2);
                        clear Aanomalia_sa
                        Pg=(GM./r.^2.*sqrt(Pg))*10^5;
                 
                    elseif volbapar(i)==20 %Gravity disturbance
                        
                        Wr_porucha=sum((AWr_porucha*covmat).*AWr_porucha,2);
                        clear AWr_porucha
                        Wr_porucha=(GM./r.^2.*sqrt(Wr_porucha));
                        
                        Wfi_porucha=sum((AWfi_porucha*covmat).*AWfi_porucha,2);
                        clear AWfi_porucha
                        Wfi_porucha=(GM./r.^2.*sqrt(Wfi_porucha));
                        
                        Wlambda_porucha=sum((AWlambda_porucha*covmat).*AWlambda_porucha,2);
                        clear AWlambda_porucha
                        Wlambda_porucha=(GM./(r.^2.*cos(fiG)).*sqrt(Wlambda_porucha));

                        Pg=sqrt(Wr_porucha.^2+Wlambda_porucha.^2+Wfi_porucha.^2)*10^5;

                        clear Wr_porucha Wfi_porucha Wlambda_porucha
                        
                    elseif volbapar(i)==21 %Gravity disturbance sa
                        
                        Pg=sum((Aporucha_sa*covmat).*Aporucha_sa,2);
                        clear Aporucha_sa
                        Pg=(GM./r.^2.*sqrt(Pg))*10^5;
                        
                    elseif volbapar(i)==22 %Height anomaly Ell
                        
                        Pg=sum((AzetaEl*covmat).*AzetaEl,2);
                        clear AzetaEl
                        Pg=(GM./(r.*gamaP).*sqrt(Pg));
     
                    elseif volbapar(i)==23 %Height anomaly
                        
                        Pg=sum((Azeta*covmat).*Azeta,2);
                        clear Azeta
                        Pg=(GM./(r.*gamaP).*sqrt(Pg));
                        
                    elseif volbapar(i)==24 %Second radial derivative of disturbing potential
                        
                        Pg=sum((AT_rr*covmat).*AT_rr,2);
                        clear AT_rr
                        Pg=(GM./r.^3.*sqrt(Pg))*10^9;
                        
                    elseif volbapar(i)==25 %Second radial derivative of gravity potential
                        
                        Pg=sum((AWrr*covmat).*AWrr,2);
                        clear AWrr
                        Pg=(GM./r.^3.*sqrt(Pg))*10^9;
                        
                    end

                    if i==1
                        P=Pg;
                        clear Pg
                    else
                        P=[P Pg];
                        if i==pocetpar
                           clear Pg
                        end
                    end                    
                end 
                
                clear q q2 gamaP covmat
            end
            
            set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow;

            cas=toc; %Stop clock

            %% Export data
            %==============================================================
                        
            if coord==1 %Spherical coordinates
                fi=fiG; %Ellipsoidal latitude is replaced by the spherical one
                clear fiG %Redundand vector fiG is deleted
            elseif coord==0 %Ellipsoidal coordinates
                clear fiG %Spherical latitude is deleted, since exported is
                %to be the ellipsoidal one
            end
            
            if get(findobj('tag','export'),'value')==1 || get(findobj('tag','datamat'),'value')==1
                    
                if volbagridcheck==1
                    if STD==0
                        [lambda,fi]=meshgrid(rad2deg(lambda),rad2deg(fi));
                        [j_fi,i_fi]=size(fi);
                        fi=fi(:);
                        lambda=lambda(:);
                    elseif STD==1
                        fi=rad2deg(fi);
                        lambda=rad2deg(lambda);
                    end
                elseif volbadiskcheck==1
                    fi=rad2deg(fi);
                    lambda=rad2deg(lambda);
                end                
                 
                set(findobj('tag','hlasky'),'string',...
                'Creating data file...','fontsize',8,...
                'foregroundcolor','k'); drawnow;
                              
                try
                    if volbagridcheck==1
                        output=[fi';lambda';P'];
                    elseif volbadiskcheck==1
                        if coord==1 %Spherical coordinates
                            output=[fi';lambda';r';P'];
                        elseif coord==0 %Ellipsoidal coordinates
                            output=[fi';lambda';h';P'];
                        end
                    end

                    if get(findobj('tag','export'),'value')==1 %Export to data file
                        [rows_output,~]=size(output);

                        exp1=fopen([outadresar,[outname '.txt']],'w');
                        fprintf(exp1,...
                            [repmat('% 0.12e ',1,rows_output-1) '% 0.12e\n'],output);                   
                        fclose(exp1);
                        set(findobj('tag','hlasky'),'string',...
                        '','fontsize',8,...
                        'foregroundcolor','k'); drawnow;  
                    end
    
                    if get(findobj('tag','datamat'),'value')==1 %Export to .mat file
                        output=output';
                        save([outadresar,[outname '.mat']],'output');
                    end
                
                    clear output
                    
                    export_error=0;
                catch err
                    export_error=1;
                    
                    save([outadresar,[outname '_phi.mat']],'fi');
                    save([outadresar,[outname '_lambda.mat']],'lambda');
                    save([outadresar,[outname '_functionals.mat']],'P');
                end
            end
            
            %% Export report
            %============================================================== 
            
            if get(findobj('tag','report'),'value')==1
                
                date_time=clock;
                date_time=datestr(date_time);
                date_time(1:11);
                date_time(13:20);
                
                reportname=[outname '_Report' '.txt'];
                exp=fopen([outadresar,reportname],'w');
                
                set(findobj('tag','hlasky'),'string',...
                'Creating report file...','fontsize',8,...
                'foregroundcolor','k'); drawnow;
                
                if get(findobj('tag','export'),'value')==0 && get(findobj('tag','datamat'),'value')==0
                    fi=rad2deg(fi);
                    lambda=rad2deg(lambda);
                    if STD==0
                        j_fi=length(fi);
                        i_fi=length(lambda);                    
                    end
                end  

                fprintf(exp,'Software                                           \tGrafLab 1.00\n');
                fprintf(exp,'Generating date                                    \t');
                fprintf(exp,'%s',date_time(1:11));
                fprintf(exp,'\nGenerating time                                  \t');
                fprintf(exp,'%s',date_time(13:20));
                fprintf(exp,'\nComputed                                         \t');
                
                if STD==0
                    fprintf(exp,'%s','Functionals of the geopotential');
                elseif STD==1
                    fprintf(exp,'%s','Commission errors');
                end
                
                if STD==0
                    fprintf(exp,'\nGeopotential model file                          \t');
                    fprintf(exp,'%s',GGMname);
                elseif STD==1
                    fprintf(exp,'\nVariance-covariance matrix file                  \t');
                    fprintf(exp,'%s',GGMcovname);
                end
                
                if STD==0
                    if get(findobj('tag','P1'),'value')==10 || get(findobj('tag','P2'),'value')==10 || get(findobj('tag','P3'),'value')==10 || get(findobj('tag','P4'),'value')==10 || get(findobj('tag','P1'),'value')==23 || get(findobj('tag','P2'),'value')==23 || get(findobj('tag','P3'),'value')==23 || get(findobj('tag','P4'),'value')==23
                        fprintf(exp,'\nDigital terain model file                          \t');
                        fprintf(exp,'%s',loadnameDMR);
                        fprintf(exp,'\nNewtonian gravitational constant (m^3.kg^-1.s^-2)  \t');
                        fprintf(exp,'%1.5e',G);
                        fprintf(exp,'\nDensity of the crust (kg.m^-3)                     \t');
                        fprintf(exp,'%4.0f',ro);
                    end
                end
                fprintf(exp,'\nGM of the geopotential model (m^3.s^-2)          \t');                
                fprintf(exp,'%1.9e',GM);
                fprintf(exp,'\nR of the geopotential model (m)                  \t');
                fprintf(exp,'%1.9e',R);
                fprintf(exp,'\nMinimum used degree                              \t');
                if nmin0==1
                    fprintf(exp,'%0.0f',0);
                elseif nmin1==1
                    fprintf(exp,'%0.0f',1);
                else
                    fprintf(exp,'%0.0f',nmin);
                end
                fprintf(exp,'\nMaximum used degree                              \t');
                fprintf(exp,'%0.0f',nmax);  
                fprintf(exp,'\nReference ellipsoid                              \t');

                if get(findobj('tag','ell'),'value') == 1;
                    fprintf(exp,'GRS80');
                elseif get(findobj('tag','ell'),'value') == 2;
                    fprintf(exp,'WGS84');
                end
                                
                fprintf(exp,'\nType of the input coordinates                    \t');
                if coord==1
                    fprintf(exp,'Spherical');
                elseif coord==0
                    fprintf(exp,'Ellipsoidal');
                end
                fprintf(exp,'\nLatitude limit North (deg)                       \t');
                fprintf(exp,'%0.9f',max(fi));
                fprintf(exp,'\nLatitude limit South (deg)                       \t');
                fprintf(exp,'%0.9f',min(fi));
                fprintf(exp,'\nLongitude limit West (deg)                       \t');
                fprintf(exp,'%0.9f',min(lambda));
                fprintf(exp,'\nLongitude limit East (deg)                       \t');
                fprintf(exp,'%0.9f',max(lambda));

                if volbagridcheck==1
                    fprintf(exp,'\nLatitude parallels                               \t');
                    fprintf(exp,'%0.0f',j_fi);
                    fprintf(exp,'\nLongitude parallels                              \t');
                    fprintf(exp,'%0.0f',i_fi);
                    fprintf(exp,'\nNumber of grid points                            \t');
                    fprintf(exp,'%0.0f',j_fi*i_fi);                                     
                    fprintf(exp,'\nGrid height above the reference surface (m)      \t');
                    
                    if coord==1 %Spherical coordinates
                        fprintf(exp,'%0.3f',hsph);
                    elseif coord==0 %Ellipsoidal coordinates
                        fprintf(exp,'%0.3f',h(1));
                    end
                    
                    fprintf(exp,'\nComputation time (s)                             \t'); 
                    fprintf(exp,'%0.0f',cas);
                    fprintf(exp,'\nComputation of fully normalized ALFs             \t');
                    
                    if volbaALFs==1
                        fprintf(exp,'Standard forward column method');
                    elseif volbaALFs==2
                        fprintf(exp,'Modified forward column method');
                    elseif volbaALFs==3
                        fprintf(exp,'Extended-range arithmetic');
                    end

                    fprintf(exp,'\n\nExported data file contains the following columns:');
                    
                    if coord==1 %Spherical coordinates
                        fprintf(exp,'\nSpherical Latitude (deg)\t');
                    elseif coord==0 %Ellipsoidal coordinates
                        fprintf(exp,'\nEllipsoidal Latitude (deg)\t');
                    end
                    fprintf(exp,'Longitude (deg)\t');                               

                elseif volbadiskcheck==1
                    fprintf(exp,'\nNumber of points                                 \t');
                    fprintf(exp,'%0.0f',length(fi));
                    fprintf(exp,'\nComputation time (s)                               \t');
                    fprintf(exp,'%0.0f',cas); 
                    fprintf(exp,'\nComputation of fully normalized ALFs             \t');
                    
                    if volbaALFs==1
                        fprintf(exp,'Standard forward column method');
                    elseif volbaALFs==2
                        fprintf(exp,'Modified forward column method');
                    elseif volbaALFs==3
                        fprintf(exp,'Extended-range arithmetic');
                    end
                    
                    fprintf(exp,'\n\nExported data file contains the following columns:');
                    
                    if coord==1 %Spherical coordinates
                        fprintf(exp,'\nSpherical Latitude (deg)\t');
                    elseif coord==0 %Ellipsoidal coordinates
                        fprintf(exp,'\nEllipsoidal Latitude (deg)\t');
                    end
                    fprintf(exp,'Longitude (deg)\t');
                    
                    if coord==1
                        fprintf(exp,'Spherical radius (m)\t');
                    elseif coord==0
                        fprintf(exp,'Height above the reference ellipsoid (m)\t');
                    end
                    
                end
                               
                for i=1:pocetpar
                    if volbapar(i)~=1
                        Par=cellstr(get(findobj('tag',sprintf('P%0.0d',i)),'string'));
                        Par=Par(volbapar(i));

                        fprintf(exp,'%s\t',char(Par));
                        
                        if volbapar(i)==2 || volbapar(i)==3 
                            fprintf(exp,'(arcsec)\t');
                        elseif volbapar(i)==4
                            fprintf(exp,'(arcsec)\tAzimuth\t(deg)\t');
                        elseif volbapar(i)==5 || volbapar(i)==11 || volbapar(i)==18
                            fprintf(exp,'(m^2.s^-2)\t');
                        elseif volbapar(i)==10 || volbapar(i)==22 || volbapar(i)==23
                            fprintf(exp,'(m)\t');
                        elseif volbapar(i)==16 || volbapar(i)==17 || volbapar(i)==19 || volbapar(i)==20 || volbapar(i)==21                        
                            fprintf(exp,'(mGal)\t');
                        elseif volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15 || volbapar(i)==24 || volbapar(i)==25
                            fprintf(exp,'(Eotvos)\t');
                        end
                    end
                end
                
                if get(findobj('tag','export'),'value')==0 && get(findobj('tag','datamat'),'value')==0
                    export_error=0;
                end
                
                if export_error==1
                    fprintf(exp,'\n\nNote: GrafLab failed to create the data file due to the lack of memory.');
                    fprintf(exp,'\nHowever, GrafLab created three *.mat files, which contain the all computed data.');
                end
                
                fclose(exp); 
                
                set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow; 
            end
            
            %% Display data
            %==============================================================
            
            if display_data==1 && volbagridcheck==1                                
                
                if STD==0
                    if get(findobj('tag','export'),'value')==1 || get(findobj('tag','datamat'),'value')==1
                        fiGrid=reshape(fi,j_fi,i_fi);
                        lambdaGrid=reshape(lambda,j_fi,i_fi);
                    elseif get(findobj('tag','report'),'value')==1
                        [lambdaGrid,fiGrid]=meshgrid(lambda,fi);
                    else
                        fi=rad2deg(fi);
                        lambda=rad2deg(lambda);
                        [lambdaGrid,fiGrid]=meshgrid(lambda,fi);
                        [j_fi,i_fi]=size(fiGrid);  
                    end
                elseif STD==1
                    if get(findobj('tag','export'),'value')==1 || get(findobj('tag','datamat'),'value')==1
                        fiGrid=reshape(fi,j_fi,i_fi);
                        lambdaGrid=reshape(lambda,j_fi,i_fi);
                    elseif get(findobj('tag','report'),'value')==1
                        fiGrid=reshape(fi,j_fi,i_fi);
                        lambdaGrid=reshape(lambda,j_fi,i_fi);
                    else
                        fi=rad2deg(fi);
                        lambda=rad2deg(lambda);
                        fiGrid=reshape(fi,j_fi,i_fi);
                        lambdaGrid=reshape(lambda,j_fi,i_fi);
                    end
                end
                
                clear fi lambda
                
                volbaformat=get(findobj('tag','nmin'),'userdata');
                colmap=get(findobj('tag','nmax'),'userdata');

                %Colormap
                if colmap==1
                    colmap='jet';
                elseif colmap==2
                    colmap='hsv';
                elseif colmap==3
                    colmap='hot';
                elseif colmap==4
                    colmap='cool';
                elseif colmap==5
                    colmap='spring';
                elseif colmap==6
                    colmap='summer';
                elseif colmap==7
                    colmap='autumn';
                elseif colmap==8
                    colmap='winter';
                elseif colmap==9
                    colmap='gray';
                elseif colmap==10
                    colmap='bone';
                elseif colmap==11
                    colmap='copper';
                elseif colmap==12
                    colmap='pink';
                elseif colmap==13
                    colmap='lines';
                end
                
                %Graphic format file
                if volbaformat==1
                    format='bmp';
                    dformat='-dbmp16m';
                elseif volbaformat==2
                    format='emf';
                    dformat='-dmeta';
                elseif volbaformat==3
                    format='eps';
                    dformat='-depsc';
                elseif volbaformat==4
                    format='jpeg';
                    dformat='-djpeg';
                elseif volbaformat==5
                    format='pdf';
                    dformat='-dpdf';
                elseif volbaformat==6
                    format='png';
                    dformat='-dpng';
                elseif volbaformat==7
                    format='tiff';
                    dformat='-dtiff';
                end
                
                coast=load('coast.mat'); %Loading continents  
                
                zobr=1;
      
                for i=1:pocetpar
                                       
                    if i==1
                        set(findobj('tag','hlasky'),'string',...
                            'Displaying 1st functional...','fontsize',...
                            8,'foregroundcolor','k'); drawnow;
                    elseif i==2
                        set(findobj('tag','hlasky'),'string',...
                            'Displaying 2nd functional...','fontsize',...
                            8,'foregroundcolor','k'); drawnow;
                    elseif i==3
                        set(findobj('tag','hlasky'),'string',...
                            'Displaying 3rd functional...','fontsize',...
                            8,'foregroundcolor','k'); drawnow;
                    elseif i==4
                        set(findobj('tag','hlasky'),'string',...
                            'Displaying 4th functional...','fontsize',...
                            8,'foregroundcolor','k'); drawnow;
                    end                    
                                                                              
                    if volbapar(i)~=1                                     
                        if volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15                   
                            for j=0:2
                                Pdisp=reshape(P(:,zobr+j),j_fi,i_fi);

                                if volbapar(i)==6
                                    nazov=['Disturbing tensor - Trr';...
                                        'Disturbing tensor - Tpp';...
                                        'Disturbing tensor - Tll'];
                                elseif volbapar(i)==7
                                    nazov=['Disturbing tensor - Trp';...
                                        'Disturbing tensor - Trl';...
                                        'Disturbing tensor - Tpl'];
                                elseif volbapar(i)==8
                                    nazov=['Disturbing tensor - Txx';...
                                        'Disturbing tensor - Tyy';...
                                        'Disturbing tensor - Tzz'];
                                elseif volbapar(i)==9
                                    nazov=['Disturbing tensor - Txy';...
                                        'Disturbing tensor - Txz';...
                                        'Disturbing tensor - Tyz'];
                                elseif volbapar(i)==12
                                    nazov=['Gravitational tensor - Vrr';...
                                        'Gravitational tensor - Vpp';...
                                        'Gravitational tensor - Vll'];
                                elseif volbapar(i)==13
                                    nazov=['Gravitational tensor - Vrp';...
                                        'Gravitational tensor - Vrl';...
                                        'Gravitational tensor - Vpl'];
                                elseif volbapar(i)==14
                                    nazov=['Gravitational tensor - Vxx';...
                                        'Gravitational tensor - Vyy';...
                                        'Gravitational tensor - Vzz'];
                                elseif volbapar(i)==15
                                    nazov=['Gravitational tensor - Vxy';...
                                        'Gravitational tensor - Vxz';...
                                        'Gravitational tensor - Vyz'];
                                end

                                Okno=figure('numbertitle','off','name',...
                                    char(nazov(j+1,:)),'visible','off');
                                worldmap([min(min(fiGrid)) max(max(fiGrid))],...
                                    [min(min(lambdaGrid)) max(max(lambdaGrid))]);
                                contourfm(fiGrid,lambdaGrid,Pdisp,ncolor,...
                                    'linestyle','none');
                                colormap(sprintf('%s(%d)',colmap,ncolor));
                                labelcolbar=colorbar('location','southoutside');

                                %Units in colorbar
                                set(get(labelcolbar,'xlabel'),'string','Eotvos');
                                
                                geoshow(coast.lat,coast.long,'color','black');
                                
                                print(Okno,sprintf('%s',dformat),...
                                    sprintf('-r%i',DPI),...
                                    sprintf('%s\\%s_%s.%s',outadresar,...
                                    outname,char(nazov(j+1,:)),format));
                            end                                                       
                            
                            zobr=zobr+3;                          
                        else                                
                            Pdisp=reshape(P(:,zobr),j_fi,i_fi);
                            
                            %Deleting azimuth if Theta has been computed
                            if volbapar(i)==4 && STD==0
                                P(:,zobr+1)=[];
                            end
                            
                            Nazov_okna=cellstr(get(findobj('tag','P1'),'string'));
                            Nazov_okna=Nazov_okna(volbapar(i));                            
                            Okno=figure('numbertitle','off','name',...
                                char(Nazov_okna),'visible','off');
                            worldmap([min(min(fiGrid)) max(max(fiGrid))],...
                                [min(min(lambdaGrid)) max(max(lambdaGrid))]);
                            contourfm(fiGrid,lambdaGrid,Pdisp,ncolor,...
                                'linestyle','none');
                            colormap(sprintf('%s(%d)',colmap,ncolor));
                            labelcolbar=colorbar('location','southoutside');
                          
                            %Units in colorbar
                            if volbapar(i)==2 || volbapar(i)==3 || volbapar(i)==4
                                set(get(labelcolbar,'xlabel'),'string','arcsec');
                            elseif volbapar(i)==5 || volbapar(i)==11 || volbapar(i)==18
                                set(get(labelcolbar,'xlabel'),'string','m^2 \cdot s^{-2}');
                            elseif volbapar(i)==10 || volbapar(i)==22 || volbapar(i)==23
                                set(get(labelcolbar,'xlabel'),'string','m');
                            elseif volbapar(i)==16 || volbapar(i)==17 || volbapar(i)==19 || volbapar(i)==20 || volbapar(i)==21                        
                                set(get(labelcolbar,'xlabel'),'string','mGal');
                            elseif volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15 || volbapar(i)==24 || volbapar(i)==25
                                set(get(labelcolbar,'xlabel'),'string','Eotvos');
                            end

                            geoshow(coast.lat,coast.long,'color','black');   
                            
                            print(Okno,sprintf('%s',dformat),...
                                sprintf('-r%i',DPI),sprintf('%s\\%s_%s.%s',...
                                outadresar,outname,char(Nazov_okna),format));
                            
                            zobr=zobr+1;
                        end  
                    end                     
                end
    
                set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow;
            end 
            
            set(findobj('tag','R_text'),'userdata','');
            set(findobj('tag','ell_text'),'userdata','');
          
            clear all

            set(findobj('tag','hlasky'),'string',...
                    'Computation has been finished','fontsize',8,...
                    'foregroundcolor','k'); drawnow; 
                
        case('Close')
            close all
    end
end