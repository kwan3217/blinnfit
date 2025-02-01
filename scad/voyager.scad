/* Notes:

Bus is not centered on rotation axis (1.12 inch offset towards +Y)
CG of spacecraft is not on rotation axis (

*/


module decagon() {
  //main body of spacecraft, with no details. 10-bay structure, numbered
  //from bay 1 perpendicular to +X. Science boom attached to bays 3 and 4,
  //Gold Record to bay 6, RTGs to bays 8 and 9, star trackers to bay 10
  translate([0,1.12,0])
    rotate([0,0,90])
    difference() {
    cylinder(h=17,r=38,$fn=10);
      translate([0,0,-1])
    cylinder(h=19,r=30,$fn=10);
  }
}

module fueltank() {
    translate([0,0,12])
    sphere(r=14);
}

decagon();
fueltank();