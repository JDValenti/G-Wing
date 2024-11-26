# -*- coding: utf-8 -*-

switchSettings = """

// NOTE: IF BOS2 not automatical copied into Wings/<<ExampleWing>>/WingOpenSCAD/, 
//       then the library can be found at: CodeFiles/OpenSCAD/BOSL2
//       copy BOSL2 Folder to Wings/<<ExampleWing>>/WingOpenSCAD/

SpanSwitch = 0; // 0 = FullSpan, 1 = SemiSpan


// Main switches, set ONE of the following to 1
WingSolidSwitch     = 1;
WingSkinSwitch      = 0;
StuctureSwitch      = 0;
Skin_w_StructSwitch = 0;

// Molds and Fillets from Old Project, not often used
USfilletSwitch      = 0;
USskinSwitch        = 0;
USskinWstructSwitch = 0;
USmoldSwitch        = 0;

LSfilletSwitch      = 0;
LSskinSwitch        = 0;
LSskinWstructSwitch = 0;
LSmoldSwitch        = 0;

// Wing OML and Internal Structure Parameters
"""

baseScript = """
w_e = t_sp;
t = w_e;
h_Struct = b/10;
n = len(WingShape);

// Mold Settings
moldHeight = 30;
moldOffset = 0;

// ---------------------------------------------
// Wing Assembly -------------------------------
// ---------------------------------------------
include <BOSL2/std.scad>
// Full Wing -------------------
if (WingSolidSwitch){ // Solid Wing
    color("grey",0.3)
    Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords);
}

if (WingSkinSwitch){ // Skin

difference(){
WingOffset(0.0*t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords,0.5*t_TE);
    scale([1,(b-t)/b,1])
    Wing(0.0*t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords);
    
}
}

if (StuctureSwitch){ // Structure Only
//translate([b/5,0,0])
    color("white")
intersection(){
curvilinearSparModule(b,xcBV,w_e,h_Struct,WingShape,WingStruct);
Wing(t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords);}
} 

if (Skin_w_StructSwitch){ // Skin w/ Structure
union(){
intersection(){
curvilinearSparModule(b,xcBV,w_e,h_Struct,WingShape,WingStruct);
Wing(t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords);}
difference(){
WingOffset(0.0*t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords,0.5*t_TE);
    scale([1,(b-t)/b,1])
    Wing(0.0*t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords);
}
}}


// Upper Surface --------------------
if (USfilletSwitch){ // US fillet Solid
    Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,USfillet);
}

if (USskinSwitch){ // US Skin
difference(){
    Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,USfillet);
    translate([0,-t,0])
    WingOffset(t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords,-t_sk);
}}

if (USskinWstructSwitch){ // US Skin w/Structure
union(){
intersection(){
translate([-xcBV*latWidth[0],0,-0.5*h_Struct])
curvilinearSparModule(b,xcBV,w_e,h_Struct,WingShape,WingStruct);
Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,USfillet);}
difference(){
    Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,USfillet);
    translate([0,-t,0])
    WingOffset(t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords,-t_sk);
}
}}

if (USmoldSwitch){ // US mold
Mold(t_TE,b,xcBV,SpanSwitch,WingShape,UScoords,moldOffset,moldHeight);
}

// Lower Surface ---------------------
if (LSfilletSwitch){ // LS fillet Solid
    Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,LSfillet);
}

if (LSskinSwitch){ // LS Skin w/o Structure
difference(){
    Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,LSfillet);
    translate([0,-t,0])
    WingOffset(t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords,-t_sk);
}}

if (LSskinWstructSwitch){ // LS Skin w/Structure
union(){
intersection(){
translate([-xcBV*WingShape[n/2][1],0,0])
curvilinearSparModule(b,xcBV,w_e,h_Struct,WingShape,WingStruct);
Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,LSfillet);}
difference(){
    Wing(0.5*t_TE,b,xcBV,SpanSwitch,WingShape,LSfillet);
    translate([200,0,0])
    translate([0,-t,0])
WingOffset(t_TE,b,xcBV,SpanSwitch,WingShape,AFcoords,-t_sk);
}
}}

if (LSmoldSwitch){ // LS mold
Mold(t_TE,b,xcBV,SpanSwitch,WingShape,LScoords,moldOffset,moldHeight);
}


// ------------------------------------------ //
// ---------------- Modules ----------------- //
// ------------------------------------------ //
module curvilinearSparModule(b,x_c_LL,w_e,h_struct,WingShape,WingStruct){

n_b  = len(WingShape);
n_sp = len(WingStruct[0]);

union(){
for(i = [0:n_sp-1]){

for(j = [0:n_b-2]){
   if(WingStruct[j][i] != 0 || WingStruct[j+1][i] != 0){
   hull(){
translate([b*WingShape[j][1]*(WingStruct[j][i]-x_c_LL) ,b*WingShape[j][0],0])
linear_extrude(h_Struct,center = true)
circle(d = w_e, $fn = 20);
translate([b*WingShape[j+1][1]*(WingStruct[j+1][i]-x_c_LL) , b*WingShape[j+1][0],0])
linear_extrude(h_Struct,center = true)
circle(d = w_e, $fn = 20);
   }
}    
}
}
}
}

module Wing(t_TE,b,x_c_BV,SpanSwitch,WingShape,SectCoords){
//include <BOSL2/std.scad>
// Measure Data Lengths ------------------------
n_f = len(SectCoords);     //#ofPoints on airfoil
n_b = len(WingShape);   //#ofPoints on wing

// Sort Wing Data ------------------------------
y   = b*[for(i = [0:n_b-1])WingShape[i][0]];
c   = b*[for(i = [0:n_b-1])WingShape[i][1]];
alpha = [for(i = [0:n_b-1])(180/3.14159)*WingShape[i][2] + alpha0L];

P_BV = [for(k = [0:n_b-1]) [0, 0, y[k]]];

surf = Surf(SectCoords);

wing = [for(j = [0:n_b-1]) fPosition(fRotate(fFiTE(c[j]*SectCoords,surf,t_TE,c[j]),c[j],alpha[j]),c[j],P_BV[j])];


rotate(a = 90,v = [1,0,0])
if (SpanSwitch == 0){
    union(){ for(i = [0:n_b-2]) loft(wing[i],wing[i+1],1);}
} else{
    difference(){
        union(){ for(i = [0:n_b/2-1]) loft(wing[i],wing[i+1],1);};
        translate(-1.2*b*[0.25,0.5,0])
        cube(1.2*b*[1,1,1]);
    }
}

function Surf(Coords) = [for (i = [0:n_f-1]) 
    (i == 0) ? 1 
    : (SectCoords[i][0]<SectCoords[i-1][0]) ? 1 
    : -1];
    // Takes in AF Coords and labels US or LS for each point

function fFiTE(Coords,surf,t_TE,c) = [for (i = [0:n_f-1]) 
    [Coords[i][0], 
     Coords[i][1] + surf[i]*0.5*(t_TE - (Coords[0][0]-Coords[n_f-1][0]))*Coords[i][0]/c]];
    // Adjust TE thickness

function fRotate(Coords,c,alpha) = [for (i = [0:n_f-1])
    Coords[i]*[[cos(alpha),-sin(alpha)],[sin(alpha),cos(alpha)]]
    + [0,0.25*c*sin(alpha)]
];

function fPosition(Coords,c,P_BV) = [for (i = [0:n_f-1])
    [Coords[i][0] - x_c_BV*c + P_BV[0],Coords[i][1] + P_BV[1],P_BV[2]]
];

module loft(upper_points, lower_points, number_of_layers)   
    polyhedron( 
    points = [
        for (i = [0 : number_of_layers])
            for (j = [0 : len(upper_points) - 1])
                [((upper_points[j][0] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][0] * i / number_of_layers)),
                ((upper_points[j][1] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][1] * i / number_of_layers)),
                ((upper_points[j][2] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][2] * i / number_of_layers))]
    ],
    faces = [
        [for (i= [0 : len(upper_points)-1]) i], // Upper plane.
        for (i = [0 : number_of_layers -1])
            for (j = [0 : len(upper_points) - 1]) // Towards lower points.
                [len(upper_points) * i + (j+1)%len(upper_points), 
                len(upper_points) * i + j, 
                len(upper_points) * (i+1) + j],
        for (i = [1 : number_of_layers])
            for (j = [0 : len(upper_points) - 1]) // Towards upper points.
                [len(upper_points) * i + j, 
                len(upper_points) * i + (j+1) % len(upper_points), 
                len(upper_points) * (i-1) + (j+1) % len(upper_points)],
        [for (i= [len(upper_points) * (number_of_layers+1) -1  : -1 : len(upper_points) * number_of_layers ]) i], // Lower plane.
    ]
);}


// ----------------------------------------------------- //
module WingOffset(t_TE,b,x_c_BV,SpanSwitch,WingShape,SectCoords,AFoffset){
//include <BOSL2/std.scad>
// Measure Data Lengths ------------------------
n_f = len(SectCoords);     //#ofPoints on airfoil
n_b = len(WingShape);   //#ofPoints on wing

// Sort Wing Data ------------------------------
y   = b*[for(i = [0:n_b-1])WingShape[i][0]];
c   = b*[for(i = [0:n_b-1])WingShape[i][1]];
alpha = [for(i = [0:n_b-1])(180/3.14159)*WingShape[i][2] + alpha0L];

P_BV = [for(k = [0:n_b-1]) [0, 0, y[k]]];

surf = Surf(SectCoords);

wing = [for(j = [0:n_b-1]) fPosition(fRotate(offset(fFiTE(c[j]*SectCoords,surf,t_TE,c[j]),delta=AFoffset, closed=true),c[j],alpha[j]),c[j],P_BV[j])];

rotate(a = 90,v = [1,0,0])
if (SpanSwitch == 0){
    union(){ for(i = [0:n_b-2]) loft(wing[i],wing[i+1],1);}
} else{
    difference(){
        union(){ for(i = [0:n_b/2-1]) loft(wing[i],wing[i+1],1);};
        translate(-1.2*b*[0.25,0.5,0])
        cube(1.2*b*[1,1,1]);
    }
}

function Surf(Coords) = [for (i = [0:n_f-1]) 
    (i == 0) ? 1 
    : (SectCoords[i][0]<SectCoords[i-1][0]) ? 1 
    : -1];
    // Takes in AF Coords and labels US or LS for each point

function fFiTE(Coords,surf,t_TE,c) = [for (i = [0:n_f-1]) 
    [Coords[i][0], 
     Coords[i][1] + surf[i]*0.5*(t_TE - (Coords[0][0]-Coords[n_f-1][0]))*Coords[i][0]/c]];
    // Adjust TE thickness

function fRotate(Coords,c,alpha) = [for (i = [0:n_f-1])
    Coords[i]*[[cos(alpha),-sin(alpha)],[sin(alpha),cos(alpha)]]
    + [0,0.25*c*sin(alpha)]
];

function fPosition(Coords,c,P_BV) = [for (i = [0:n_f-1])
    [Coords[i][0] - x_c_BV*c + P_BV[0],Coords[i][1] + P_BV[1],P_BV[2]]
];

module loft(upper_points, lower_points, number_of_layers)   
    polyhedron( 
    points = [
        for (i = [0 : number_of_layers])
            for (j = [0 : len(upper_points) - 1])
                [((upper_points[j][0] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][0] * i / number_of_layers)),
                ((upper_points[j][1] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][1] * i / number_of_layers)),
                ((upper_points[j][2] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][2] * i / number_of_layers))]
    ],
    faces = [
        [for (i= [0 : len(upper_points)-1]) i], // Upper plane.
        for (i = [0 : number_of_layers -1])
            for (j = [0 : len(upper_points) - 1]) // Towards lower points.
                [len(upper_points) * i + (j+1)%len(upper_points), 
                len(upper_points) * i + j, 
                len(upper_points) * (i+1) + j],
        for (i = [1 : number_of_layers])
            for (j = [0 : len(upper_points) - 1]) // Towards upper points.
                [len(upper_points) * i + j, 
                len(upper_points) * i + (j+1) % len(upper_points), 
                len(upper_points) * (i-1) + (j+1) % len(upper_points)],
        [for (i= [len(upper_points) * (number_of_layers+1) -1  : -1 : len(upper_points) * number_of_layers ]) i], // Lower plane.
    ]
);}

// ----------------------------------------------------- //
module Mold(t_TE,b,x_c_BV,SpanSwitch,WingShape,SurfCoords,AFoffset,moldHeight){

// Measure Data Lengths ------------------------
n_f = len(SurfCoords);     //#ofPoints on airfoil
n_b = len(WingShape);   //#ofPoints on wing

surfSwitch = (SurfCoords[0][1]<=SurfCoords[1][1])? 1 : -1;
// Sort Wing Data ------------------------------
y   = b*[for(i = [0:n_b-1])WingShape[i][0]];
c   = b*[for(i = [0:n_b-1])WingShape[i][1]];
alpha = [for(i = [0:n_b-1])(180/3.14159)*WingShape[i][2] + alpha0L];

P_BV = [for(k = [0:n_b-1]) [0, 0, y[k]]];

wingFiTE   = [for(j = [0:n_b-1]) fFiTE(c[j]*SurfCoords,surfSwitch,t_TE,c[j])];
wingOffset = [for(j = [0:n_b-1]) offset(wingFiTE[j],delta=AFoffset, closed=true)];
wingRotate = [for(j = [0:n_b-1]) fRotate(wingOffset[j],c[j],alpha[j])];
moldSects  = [for(j = [0:n_b-1]) concat(wingRotate[j],[[150,surfSwitch*moldHeight],[-25,surfSwitch*moldHeight]])];
wing       = [for(j = [0:n_b-1]) fPosition(moldSects[j],c[j],P_BV[j])];

rotate(a = 90,v = [1,0,0])
if (SpanSwitch == 0){
    union(){ for(i = [0:n_b-2]) loft(wing[i],wing[i+1],1);}
} else{
    difference(){
        union(){ for(i = [0:n_b/2-1]) loft(wing[i],wing[i+1],1);};
        translate(-2*max(c)*[0.5,0.5,0])
        cube(2*max(c)*[1,1,1]);
    }
}

function fFiTE(Coords,surfSwitch,t_TE,c) = [for (i = [0:n_f-1]) 
    [Coords[i][0], 
     Coords[i][1] + surfSwitch*(0.5*t_TE -Coords[n_f-1][0]/c)*Coords[i][0]/c]];
    // Adjust TE thickness

function fRotate(Coords,c,alpha) = [for (i = [0:n_f-1])
    Coords[i]*[[cos(alpha),-sin(alpha)],[sin(alpha),cos(alpha)]]
    + [0,x_c_BV*c*sin(alpha)]
];

function fPosition(Coords,c,P_BV) = [for (i = [0:n_f+1])
    [Coords[i][0] - x_c_BV*c + P_BV[0],Coords[i][1] + P_BV[1],P_BV[2]]
];

module loft(upper_points, lower_points, number_of_layers)   
    polyhedron( 
    points = [
        for (i = [0 : number_of_layers])
            for (j = [0 : len(upper_points) - 1])
                [((upper_points[j][0] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][0] * i / number_of_layers)),
                ((upper_points[j][1] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][1] * i / number_of_layers)),
                ((upper_points[j][2] * (number_of_layers - i) / number_of_layers)
                + (lower_points[j][2] * i / number_of_layers))]
    ],
    faces = [
        [for (i= [0 : len(upper_points)-1]) i], // Upper plane.
        for (i = [0 : number_of_layers -1])
            for (j = [0 : len(upper_points) - 1]) // Towards lower points.
                [len(upper_points) * i + (j+1)%len(upper_points), 
                len(upper_points) * i + j, 
                len(upper_points) * (i+1) + j],
        for (i = [1 : number_of_layers])
            for (j = [0 : len(upper_points) - 1]) // Towards upper points.
                [len(upper_points) * i + j, 
                len(upper_points) * i + (j+1) % len(upper_points), 
                len(upper_points) * (i-1) + (j+1) % len(upper_points)],
        [for (i= [len(upper_points) * (number_of_layers+1) -1  : -1 : len(upper_points) * number_of_layers ]) i], // Lower plane.
    ]
);}
$vpd = 2*b;
$vpr = [0,0,90];
"""